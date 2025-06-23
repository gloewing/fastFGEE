# Internal estimation function
#
fun.gee1step.dist <- function(orig.data, dx, formula, family, X_, Y_, namesd,
                              N_clusters, clusters, glmfit, yindex.vec, bs,
                              cov.type = "exchangeable", parallel = FALSE,
                              sandwich = FALSE, cv = TRUE, cv.grid, B = 2000,
                              rho.smooth = TRUE, joint.CI = TRUE, weighted = TRUE,
                              return.early = FALSE, V.inv = NULL, time = NULL,
                              index = "yindex.vec") {


  p <- NULL
  resid <- NULL
  wt_ij <- NULL
  rho_ij <- NULL
  sum_r <- NULL
  uss_r <- NULL
  N <- NULL
  cname_ <- "cname_" #NULL # GCL changed 3-2-25
  v <- NULL
  residv <- NULL
  Y <- NULL
  yindex.vec <- NULL
  rho <- NULL
  time <- NULL

  # if(parallel)  num_cores <- as.integer(round(parallel::detectCores() * 0.85))
  ###

  dr <- data.table::copy(dx) # for robust se

  dx[, p := stats::predict.glm(glmfit, type = "response")]

  if (family == "binomial") {
    dx[, v := p * (1-p)]
  }else if (family == "poisson") {
    dx[, v := p]
  }else if (family == "gaussian") {
    dx[, resid := as.numeric(glmfit$residuals)  ] # vectorize outcome from original data]
    dx[, v := stats::var(resid), by = yindex.vec] # this assumes homoscedasticity
    dx[, resid := NULL]
  }

  dx[, resid := (Y - p) / sqrt(v) ]
  namesd <- X_

  # correlation parameter
  if(weighted){
    rho <- corr.estimate(dx, corr = cov.type,
                         rho.smooth = rho.smooth,
                         glmfit = glmfit,
                         index = index)
    dx[, index.vec := get(index)]
    dx <- dx[rho, on = .(index.vec)] # add rho
    rho_vec <- as.numeric(rho$rho[order(rho$index.vec)])
  }else{
    rho_vec <- rep(0, length(rho_vec)) # induces working independence structure
    dx[,rho := 0]
  }

  # get penalties
  sp_vec <- glmfit["sp"][[1]]
  penalty_diag <- penalty.construct(model=glmfit, nonSmooths = 0) # nonSmooths set to 0 here
  if(all(sp_vec==0)){
    penalty_diag <- matrix(0, ncol = ncol(penalty_diag), nrow = nrow(penalty_diag))
  }

  ### Estimate beta
  yindex <- unique(dx$yindex.vec)
  L <- length(yindex)
  message("Calculating W1")
  if(is.null(V.inv)){
    wi = W.estimate(dx = dx, namesd = namesd,
                    cname_ = cname_, family = family,
                    corr = cov.type, dt = TRUE,
                    index = index)
  }else{
    # exact V.inv given
    message("full W")
    dx_list = lapply(clusters, function(i) dx[cname_ == i])
    wi <- lapply(clusters, function(i)
      .getW_true(dx_list[[i]],
                 namesd,
                 V.inv = V.inv,
                 rho = NULL,
                 family = family))
  }

  gc()
  message("Calculating D1")
  if(is.null(V.inv)){
    di <- D.estimate(dx = dx, namesd = namesd,
                     cname_ = cname_, family = family,
                     corr = cov.type, dt = TRUE,
                     index = index)
  }else{
    message("full D")
    di <- lapply(clusters, function(i)
      .getD_true(dx_list[[i]],
                 namesd,
                 V.inv = V.inv,
                 rho = NULL,
                 family = family))
    rm(dx_list)
  }
  gc()

  if(return.early){
    result <- list(beta = as.vector(glmfit$coefficients),
                   rho = rho_vec,
                   di = di,
                   wi = wi,
                   model = glmfit)

    return(result)
  }

  ##############
  # tuning
  ##############
  message("Smoothing Parameter Tuning")
  if(!is.null(cv.grid)){
    if(is.list(cv.grid)){
      grid <- cv.grid
      if(length(cv.grid) > 1){
        tune.ind <- TRUE
        # this only needs to be applied to the first one
      }else{
        tune.ind <- ifelse(length(grid[[1]]) > 1, TRUE, FALSE) # don't tune if there's only 1 tuning parameter value in the one set
        if(tune.ind){
          # if its a list but only one element (one grid) then duplicate for pre-tune
          grid <- list(grid, grid) # dublicate for pre-tune
        }
      }
    }else{
      # not a list
      grid <- list(cv.grid)
      tune.ind <- ifelse(length(grid[[1]]) > 1, TRUE, FALSE) # don't tune if there's only 1 tuning parameter set
    }

    # pre-tuning (scale smoothing parameters by same scalar - not expand.grid)
    if(tune.ind)  grid[[1]] <- do.call(cbind, lapply(glmfit$sp, function(x) x*grid[[1]]))

    if(tune.ind){
      # if multiple tuning values provided
      cv.lambda <- fun.gee1step.cv(w = wi,
                                   d = di,
                                   namesd = namesd,
                                   grid,
                                   data = dx,#orig.data,
                                   cname_ = "cluster", # may need to update this to be more general
                                   fit.initial = glmfit,
                                   cv = cv,
                                   B = 500,
                                   seed = 1)

      message("CV Complete")
      penalty_diag <- penalty.construct(model=glmfit,
                                        lambda = cv.lambda$lambda.star,
                                        nonSmooths = 0) # nonSmooths set to 0 here
    }else{
      cv.lambda <- data.frame(mse = NA,
                              lambda.star = grid,
                              grid = grid,
                              bias = NA,
                              asy.var = NA)

      penalty_diag <- penalty.construct(model=glmfit,
                                        lambda = grid,
                                        nonSmooths = 0) # nonSmooths set to 0 here
    }

  }else{
    # do no tune
    message("No penalty")
    penalty_diag <- penalty.construct(model=glmfit,
                                      lambda = glmfit$sp * 0,
                                      nonSmooths = 0)
    cv.lambda <- NULL
  }

  # bootstrap
  if(sandwich == "fastBoot" | sandwich == "boot"){
    # estimate var/cov matrix if bootstrap
    # need to use initial di, wi, beta (not after one-step)
    # IDs needed to adjust penalty for varying cluster sizes
    clust.vec <- as.vector(dx[,..cname_])[[1]]
    vb <- var.est(di = di,
                  wi = wi,
                  beta = as.numeric(glmfit$coefficients),
                  penalty_diag = penalty_diag,
                  clust.vec = clust.vec,
                  B = B,
                  sandwich = sandwich,
                  exact = FALSE)
  }

  beta <- glmfit$coefficients
  penalty_vec <- as.numeric(penalty_diag %*% beta)
  wi <- Reduce("+", wi)/N_clusters
  di <- Reduce("+", di)
  inf_fn <- sanic::solve_chol(a = wi + penalty_diag,
                              b = di - penalty_vec*N_clusters)
  beta2 <- glmfit$coefficients <- beta + inf_fn / N_clusters
  penalty_vec <- as.numeric(penalty_diag %*% beta2)
  rm(di, wi, dx)

  # --------------------------------------------
  # update wi and di using one-step estimates
  # --------------------------------------------
  ### Robust se
  dvars <- as.matrix(dr[, X_, with = FALSE])
  dr[, p := dvars %*% beta2]

  if (family == "binomial") {
    dr[, p := 1/(1 + exp(-p))]
    dr[, v := p * (1 - p)]
  }else if (family == "poisson") {
    dr[, p := exp(p)]
    dr[, v := p]
  }else if (family == "gaussian") {
    dr[, v:= stats::var( Y - p )]
  }

  dr[, resid := ( Y - p) / sqrt( v )]
  dr[, residv := ( Y - p) / v ]
  namesd <- X_

  ### Re-Estimate correlation parameter
  if(weighted){
    rho <- corr.estimate(dr, corr = cov.type,
                         rho.smooth = rho.smooth,
                         glmfit = glmfit,
                         index = index)
    dr[, index.vec := get(index)]
    dr <- dr[rho, on = .(index.vec)] # add rho
    rho_vec <- as.numeric(rho$rho)
  }else{
    rho_vec <- rep(0, length(rho_vec)) # induces working independence structure
    dr[,rho := 0]
  }

  # recalculate D and W
  message("Calculating W2")
  if(is.null(V.inv)){
    wi = W.estimate(dx = dr, namesd = namesd,
                    cname_ = cname_, family = family,
                    corr = cov.type, dt = TRUE,
                    index = index)
  }else{
    # exact V.inv given
    dr_list = lapply(clusters, function(i) dr[cname_ == i])
    wi <- lapply(clusters, function(i)
      .getW_true(dr_list[[i]],
                 namesd,
                 V.inv = V.inv,
                 rho = NULL,
                 family = family))
  }
  gc()
  message("Calculating D2")

  if(is.null(V.inv)){
    di <- D.estimate(dx = dr, namesd = namesd,
                     cname_ = cname_, family = family,
                     corr = cov.type, dt = TRUE,
                     index = index)
  }else{
    di <- lapply(clusters, function(i)
      .getD_true(dr_list[[i]],
                 namesd,
                 V.inv = V.inv,
                 rho = NULL,
                 family = family))
    rm(dr_list)
  }
  gc()

  # estimate var/cov matrix
  if(sandwich == "fastBoot" | sandwich == "boot"){ #
    message("fastBoot update")
    glmfit$Vp <- vb
  }else{
    glmfit$Vp <- vb <- var.est(di = di,
                               wi = wi,
                               beta = beta2,
                               penalty_diag = penalty_diag,
                               B = B,
                               sandwich = sandwich,
                               exact = FALSE)
  }

  ##############
  # joint CIs
  ##############
  # create matrix to generate smooth terms
  cov.names <- colnames(glmfit$model)
  y.nm <- glmfit$formula[[2]] # name of outcome in formula
  cov.names <- cov.names[!cov.names %in% c(y.nm, "yindex.vec")]
  terms.mat <- as.data.frame(matrix(1, ncol = length(cov.names), nrow = 1))
  colnames(terms.mat) <- cov.names

  if(joint.CI != FALSE){
    # Obtain qn to construct joint CI
    if(joint.CI != "np"){
      message("joint CIs: parametric bootstrap")
      qn <- joint.qn(glmfit = glmfit, NN.sim = 10000)
    }else{
      message("joint CIs: non-parametric bootstrap")
      qn <- joint.basis.np(glmfit = glmfit,
                     di = di,
                     wi = wi,
                     beta = beta2,
                     penalty_diag = penalty_diag,
                     bs = bs,
                     var.mat = glmfit$Vp,
                     exact = FALSE,
                     grid.size = 500, # num. of points to eval function on
                     NN.sim = 10000)
    }
  }

  ### Get results
  result <- list(beta = as.vector(beta2),
                 vb = vb,
                 rho = rho_vec,
                 di = di,
                 wi = wi,
                 model = glmfit,
                 pen.mat = penalty_diag,
                 lambda = cv.lambda$lambda.star,
                 grid = data.frame(mse = cv.lambda$mse,
                                   grid = cv.lambda$grid,
                                   bias = cv.lambda$bias,
                                   asy.var = cv.lambda$asy.var),
                 qn = qn)

  return(result)
}

# Internal estimation function - exact for linear outcomes
#
fun.gee1step.exact <- function(orig.data, dx, formula, family, X_, Y_, namesd, N_clusters, bs,
                              clusters, glmfit, yindex.vec, cov.type = "exchangeable", parallel = FALSE, sandwich = FALSE,
                              cv = TRUE, cv.grid, B = 2000, weighted = TRUE,
                              rho.smooth = FALSE, joint.CI = TRUE, V.inv = NULL,
                              time = NULL, index = "yindex.vec") {


  p <- NULL
  resid <- NULL
  wt_ij <- NULL
  rho_ij <- NULL
  sum_r <- NULL
  uss_r <- NULL
  N <- NULL
  cname_ <- "cname_"
  residv <- NULL
  Y <- NULL
  yindex.vec <- NULL
  rho <- NULL
  time <- NULL

  # if(parallel)  num_cores <- as.integer(round(parallel::detectCores() * 0.85))
  ###

  dr <- data.table::copy(dx) # for robust se
  dx[, p := stats::predict.glm(glmfit, type = "response")]
  dx[, resid := as.numeric(t(resid(glmfit)))  ] # vectorize outcome from original data]
  dx[, v := stats::var(resid), by = yindex.vec] # this assumes homoscedasticity
  dx[, resid := NULL]
  dx[, resid := (Y - p) / sqrt(v) ]
  namesd <- X_  # just recycle names to avoid copying

  # correlation parameter
  rho <- corr.estimate(dx, corr = cov.type,
                       rho.smooth = rho.smooth,
                       glmfit = glmfit,
                       index = index)
  dx[, index.vec := get(index)]
  dx <- dx[rho, on = .(index.vec)] # add rho
  rho_vec <- as.numeric(rho$rho[order(rho$index.vec)])

  if(!weighted){
    rho_vec <- rep(0, length(rho_vec)) # induces working independence structure
    dx[,rho := 0]
  }
  ##########################
  # consider smoothing rho instead of rho <- drho[, ...] line above (i.e., pointwise averaging)
  ##########################
  # get penalties
  sp_vec <- glmfit["sp"][[1]]
  penalty_diag <- penalty.construct(model=glmfit, nonSmooths = 0) # nonSmooths set to 0 here
  if(all(sp_vec==0)){
    penalty_diag <- matrix(0, ncol = ncol(penalty_diag), nrow = nrow(penalty_diag))
  }

  ### Estimate beta
  yindex <- unique(dx$yindex.vec)
  L <- length(yindex)
  message("Calculating W1")
  if(is.null(V.inv)){
    wi = W.estimate(dx = dx, namesd = namesd,
                    cname_ = cname_, family = family,
                    corr = cov.type, dt = TRUE,
                    index = index)
  }else{
    # exact V.inv given
      dx_list = lapply(clusters, function(i) dx[cname_ == i])
      wi <- lapply(clusters, function(i)
                 .getW_true(dx_list[[i]],
                       namesd,
                       V.inv = V.inv,
                       rho = NULL,
                       family = family))
  }

  message("Calculating D1")
  if(is.null(V.inv)){
    di <- D.estimate(dx = dx, namesd = namesd,
                           cname_ = cname_, family = family,
                           corr = cov.type, dt = TRUE,
                          exact = TRUE, index = index)
  }else{
    di <- lapply(clusters, function(i)
      .getD_true_exact(dx_list[[i]],
                 namesd,
                 V.inv = V.inv,
                 rho = NULL,
                 family = family))
    rm(dx_list)
  }

  #-----------------------------------------
  ##############
  # tuning
  ##############
  message("Smoothing Parameter Tuning")
  if(!is.null(cv.grid)){
    if(is.list(cv.grid)){
      grid <- cv.grid
      if(length(cv.grid) > 1){
        tune.ind <- TRUE
        # this only needs to be applied to the first one
      }else{
        tune.ind <- ifelse(length(grid[[1]]) > 1, TRUE, FALSE) # don't tune if there's only 1 tuning parameter value in the one set
        if(tune.ind){
          # if its a list but only one element (one grid) then duplicate for pre-tune
          grid <- list(grid, grid) # dublicate for pre-tune
        }
      }
    }else{
      # not a list
      grid <- list(cv.grid)
      tune.ind <- ifelse(length(grid[[1]]) > 1, TRUE, FALSE) # don't tune if there's only 1 tuning parameter set
    }

    # pre-tuning (scale smoothing parameters by same scalar - not expand.grid)
    if(tune.ind)  grid[[1]] <- do.call(cbind, lapply(glmfit$sp, function(x) x*grid[[1]]))

    if(tune.ind){
      # if multiple tuning values provided
      cv.lambda <- fun.gee1step.cv.exact(w = wi,
                                         d = di,
                                         namesd = namesd,
                                         grid,
                                         data = dx,#orig.data,
                                         cname_ = "cluster", # may need to update this to be more general
                                         fit.initial = glmfit,
                                         cv = cv,
                                         B = 500,
                                         seed = 1)

      message("CV Complete")
      penalty_diag <- penalty.construct(model=glmfit,
                                        lambda = cv.lambda$lambda.star,
                                        nonSmooths = 0) # nonSmooths set to 0 here
    }else{
      cv.lambda <- data.frame(mse = NA,
                              lambda.star = grid,
                              grid = grid,
                              bias = NA,
                              asy.var = NA)

      penalty_diag <- penalty.construct(model=glmfit,
                                        lambda = grid,
                                        nonSmooths = 0) # nonSmooths set to 0 here
    }
  }else{
    # do no tune
    message("No penalty")
    penalty_diag <- penalty.construct(model=glmfit,
                                      lambda = glmfit$sp * 0,
                                      nonSmooths = 0)
    cv.lambda <- NULL
  }

  # bootstrap
  if(sandwich %in% c("fastBoot", "boot") ){
    # estimate var/cov matrix if bootstrap
    # need to use initial di, wi, beta (not after one-step)
    # IDs needed to adjust penalty for varying cluster sizes
    clust.vec <- as.vector(dx[,..cname_])[[1]]
    vb <- var.est(di = di,
                  wi = wi,
                  beta = as.numeric(glmfit$coefficients),
                  penalty_diag = penalty_diag,
                  clust.vec = clust.vec,
                  B = B,
                  sandwich = sandwich,
                  exact = TRUE)
  }


  beta2 <- sanic::solve_chol(a = Reduce("+", wi) + penalty_diag,
                             b = Reduce("+", di))
  glmfit$coefficients <- beta2
  rm(di, wi, dx)

  ### Robust se
  namesd <- X_  # just recycle names to avoid copying
  dvars <- as.matrix(dr[, X_, with = FALSE])
  dr[, p := dvars %*% beta2]
  dr[, v:= stats::var( Y - p )] # gaussian
  dr[, resid := ( Y - p) / sqrt( v )]
  dr[, residv := ( Y - p) / v ]

  ### Estimate correlation parameter
  rho <- corr.estimate(dr, corr = cov.type,
                       rho.smooth = rho.smooth,
                       glmfit = glmfit,
                       index = index)
  dr[, index.vec := get(index)]
  dr <- dr[rho, on = .(index.vec)] # add rho
  rho_vec <- as.numeric(rho$rho)
  if(!weighted){
    rho_vec <- rep(0, length(rho_vec)) # induces working independence structure
    dr[,rho := 0]
  }

  # recalculate D and W
  message("Calculating W2")
  if(is.null(V.inv)){
    wi = W.estimate(dx = dr, namesd = namesd,
                    cname_ = cname_, family = family,
                    corr = cov.type, dt = TRUE)
  }else{
    # exact V.inv given
    dr_list = lapply(clusters, function(i) dr[cname_ == i])
    wi <- lapply(clusters, function(i)
                  .getW_true(dr_list[[i]],
                       namesd,
                       V.inv = V.inv,
                       rho = NULL,
                       family = family))
  }

  message("Calculating D2")
  if(is.null(V.inv)){
    di <- D.estimate(dx = dr, namesd = namesd,
                     cname_ = cname_, family = family,
                     corr = cov.type,
                     dt = TRUE,
                     exact = FALSE # FALSE b/c sandwich estimator for exact uses same di as One-Step
                     )
  }else{
    message("full D")
    di <- lapply(clusters, function(i)
      .getD_true(dr_list[[i]], # Needs this version NOT .getD_true_exact b/c sandwich uses this version
                 namesd,
                 V.inv = V.inv,
                 rho = NULL,
                 family = family))
    rm(dr_list)
  }

  # variance of coefficients
  # estimate var/cov matrix
  if(sandwich == "fastBoot" | sandwich == "boot"){ #
    message("fastBoot update")
    glmfit$Vp <- vb
  }else{
    glmfit$Vp <- vb <- var.est(di = di,
                               wi = wi,
                               beta = beta2,
                               penalty_diag = penalty_diag,
                               B = B,
                               sandwich = sandwich,
                               exact = TRUE)
  }

  if(joint.CI != FALSE){

    # Obtain qn to construct joint CI
    if(joint.CI != "np"){
      message("joint CIs: parametric bootstrap")
      qn <- joint.qn(glmfit = glmfit, NN.sim = 10000)
    }else{
      message("joint CIs: non-parametric bootstrap")
      qn <- joint.basis.np(glmfit = glmfit,
                           di = di,
                           wi = wi,
                           beta = beta2,
                           penalty_diag = penalty_diag,
                           bs = bs,
                           var.mat = glmfit$Vp,
                           exact = TRUE,
                           grid.size = 500, # num. of points to eval function on
                           NN.sim = 10000)
    }
  }else{
    qn = NULL
  }

  ### Get results
  result <- list(beta = as.vector(beta2),
                 vb = vb,
                 rho = rho_vec,
                 di = di,
                 wi = wi,
                 model = glmfit,
                 pen_mat = penalty_diag,
                 lambda = cv.lambda$lambda.star,
                 qn = qn,
                 grid = data.frame(mse = cv.lambda$mse,
                                   grid = cv.lambda$grid,
                                   bias = cv.lambda$bias,
                                   asy.var = cv.lambda$asy.var))

  return(result)
}


# sandwich
fun.gee1step.sandwich  <- function(orig.data, dx, formula, family, X_, Y_, namesd,
                                   N_clusters, clusters, glmfit, yindex.vec,
                                   cov.type = "exchangeable", parallel = FALSE,
                                   sandwich = FALSE, cv = TRUE, cv.grid, B = 2000,
                                   rho.smooth = TRUE, joint.CI = TRUE,
                                   V.inv = NULL, time = NULL, index = "yindex.vec") {


  p <- NULL
  resid <- NULL
  wt_ij <- NULL
  rho_ij <- NULL
  sum_r <- NULL
  uss_r <- NULL
  N <- NULL
  cname_ <- "cname_"
  v <- NULL
  residv <- NULL
  Y <- NULL
  yindex.vec <- NULL
  rho <- NULL
  time <- NULL

  # if(parallel)  num_cores <- as.integer(round(parallel::detectCores() * 0.85))
  ###

  dr <- data.table::copy(dx) # for robust se

  dx[, p := stats::predict.glm(glmfit, type = "response")]

  if (family == "binomial") {
    dx[, v := p * (1-p)]
  }else if (family == "poisson") {
    dx[, v := p]
  }else if (family == "gaussian") {
    dx[, resid := as.numeric(t(resid(glmfit)))  ] # vectorize outcome from original data]
    dx[, v := stats::var(resid), by = yindex.vec] # this assumes homoscedasticity
    dx[, resid := NULL]
  }

  dx[, resid := (Y - p) / sqrt(v) ]
  namesd <- X_

  # correlation parameter
  dx[, index.vec := get(index)]
  dx <- dx[,rho := 0] # fake rho
  rho_vec <- rep(0, L)
  print(rho_vec)

  # get penalties
  sp_vec <- glmfit["sp"][[1]]
  penalty_diag <- penalty.construct(model=glmfit, nonSmooths = 0) # nonSmooths set to 0 here

  ### Estimate beta
  yindex <- unique(dx$yindex.vec)
  L <- length(yindex)
  message("Calculating W1")

  wi = W.estimate(dx = dx, namesd = namesd,
                  cname_ = cname_, family = family,
                  corr = cov.type, dt = TRUE,
                  index = index)

  gc()
  message("Calculating D1")
  di <- D.estimate(dx = dx, namesd = namesd,
                   cname_ = cname_, family = family,
                   corr = cov.type, dt = TRUE,
                   index = index)


  gc()

  beta <- glmfit$coefficients
  glmfit$Vp <- vb <- var.est(di = di,
                             wi = wi,
                             beta = beta,
                             penalty_diag = penalty_diag,
                             B = B,
                             sandwich = FALSE,
                             exact = FALSE)


  ##############
  # joint CIs
  ##############
  # create matrix to generate smooth terms
  cov.names <- colnames(glmfit$model)
  y.nm <- glmfit$formula[[2]] # name of outcome in formula
  cov.names <- cov.names[!cov.names %in% c(y.nm, "yindex.vec")]
  terms.mat <- as.data.frame(matrix(1, ncol = length(cov.names), nrow = 1))
  colnames(terms.mat) <- cov.names

  if(joint.CI){
    # Obtain qn to construct joint CI
    qn <- joint.qn(glmfit = glmfit, NN.sim = 10000)
  }

  ### Get results
  result <- list(beta = as.numeric(beta),
                 vb = vb,
                 rho = rho_vec,
                 di = di,
                 wi = wi,
                 model = glmfit,
                 pen.mat = penalty_diag,
                 lambda = sp_vec,
                 grid = NULL,
                 qn = qn)

  return(result)
}

#' Estimate parameters using one-step algorithm
#' @importFrom sanic solve_chol
#' @importFrom Rfast rowMaxs ar1 colMaxs
#' @import data.table mgcv refund Matrix SuperGauss
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data a required data frame or data.table containing the variables in the model.
#' @param cluster the name of the field that identifies the clusters.
#' @param family the distribution family: gaussian, binomial, and poisson
#' @param cov.type working correlation can be "exchangeable", "ar1", or "independence. "ar1" only supported for even grids (contact package author if you need an implementation for uneven grids!)
#' @param time variable name for the time variable (necessary only for "ar1").
#' @param long.dir Logical. Defaults to \code{TRUE}, which models correlation in the longitudinal direction with the specified correlation above (cov.type). \code{TRUE} models correlation along the functional domain with an ar1 structure
#' @param sandwich Specifies the variance estimator. Set to "fastBoot" for fast cluster bootstrapping, or "boot" for standard cluster bootstrapping.
#' @param pffr.mod Instead of specifying the formula above, one can instead fit a refund::pffr() model and feed in the model object to provide the initial fit.
#' @param knots Number of knots used for initial pffr fit. Default to round(L/4) where L is the number of functional domain points.
#' @param bs Spline basis type used.
#' @param cv Cross-validation type. Defaults to cluster "fastkfold" but can also use "kfold" for standard cluster CV, or TRUE for hold-one-cluster-out CV.
#' @param cv.grid A list of vectors, where each vector is a grid to sequentially tune on. Default uses a list of three lists to tune over, which seems to work well in practice.
#' @param exact Defaults to FALSE. For continuous outcomes, setting `exact=TRUE` uses the closed-form expression (a penalized weighted-least squares) instead of the closely-related one-step estimator.
#' @param rhoSmooth Defaults to FALSE Determines wheter the pointwise correlation parameter estimates are smoothed across the functional domain.
#' @param joint.CI Defaults to TRUE which uses a parameteric bootstrap. For a non-parametric bootstrap, set `joint.CI = "np"`.
#' @param V.inv Defaults to NULL. This is an optional inverse-covariance matrix (of dimension Ln_i x Ln_i) to be used instead of the working correlations above. Experimental feature that assumes fixed sample size across subjects. Contact package authors if you need this feature.
#' @param return.early Defaults to FALSE. Calculates one-step components and returns them before parameter tuning or one-step update.
#' @param ... currently disregarded
#' @references Loewinger, G., Levis, A., Cui, E., Pereira, F. (2025). Fast Functional Generalized Estimating Equations for
#' Longitudinal Functional Data at Scale.
#' @return A "fgee1step" object. \code{model} is the updated one-step model object, so mgcv::plot.gam() can be applied to this object. \code{pffr_initial.fit} contains the initial fit refund::pffr object.
#' \code{vb} contains the variance/covariance matrix for the coefficient estimates (sandwich or bootstrap-based). \code{qn} contains the joint CI quantiles. \code{di} and \code{wi} are lists of length N, with the updated cluster-specific estimating equation and hessian terms (without the penalty).
#' @author Gabriel Loewinger \email{gloewinger@@gmail.com}
#' @examples
#' library(fastFGEE)
#' data("DTI", package = "refund")
#' set.seed(1)
#' DTI_use <- DTI[DTI$ID %in% sample(DTI$ID, 45),]
#' DTI_use <- data.frame(cca = I(DTI_use$cca),
#'                       case = DTI_use$case, visit = 1*(DTI_use$visit),
#'                       sex = DTI_use$sex,
#'                       ID = DTI_use$ID)
#'  fit_dti <- fastFGEE::fgee(formula = cca ~ case + visit,
#'                            cluster = "ID",
#'                            data = DTI_use,
#'                            family = "gaussian",
#'                            cov.type = "exchangeable")
#'  fastFGEE::fgee.plot(fit_dti)
#'  plot(fit_dti$model)
#' @export

fgee <- function(formula, data, cluster, family,
                         cov.type = c("exchangeable", "ar1", "independence"),
                         time = NULL,
                         long.dir = TRUE,
                         sandwich = TRUE,
                         pffr.mod = NULL,
                         knots = NULL,
                         bs = "bs",
                         cv = "fastkfold", cv.grid = NULL,
                         exact = FALSE, rho.smooth = FALSE,
                         joint.CI = TRUE, return.early = FALSE,
                         V.inv = NULL, ...) {

  # "declare" vars to avoid global NOTE
  orig.formula <- formula
  cname_ <- NULL
  Y <- NULL
  N <- NULL

  ### Check arguments

  if ( ! inherits(formula, "formula" ) ) {
    stop("The argument `formula` is not properly specified")
  }

  if  (! is.data.frame(data) ) {
    stop("Data must be a data.frame or data.table")
  }


  if ( ! (cluster %in% names(data)) ) {
    stop(paste("Cluster variable", cluster, "is not in data set"))
  }

  Y_nm <- all.vars(formula)[1]
  L <- ncol(data[, Y_nm])

  # fit initial pffr model
  index <- ifelse(long.dir, "yindex.vec", "time")
  if(is.null(knots))  knots <- round(L/4)
  message("Fitting pffr() model")
  # need to order by cluster, yindex.vec, time in the appropriate order
  if(is.null(pffr.mod)){
    fit_pffr <- refund::pffr(formula = formula,
                             algorithm = "bam",
                             family = family,
                             discrete = TRUE,
                             bs.yindex = list(bs = bs,
                                              k = knots,
                                              m = c(2, 1)),
                             data = data)
  }else{
    fit_pffr <- pffr.mod
  }

  yindex.vec <- fit_pffr$model$yindex.vec # functional domain index
  MM <- suppressWarnings(stats::model.matrix(fit_pffr))

  Y <- fit_pffr$model[,Y_nm]
  if(ncol(data[,Y_nm]) > 1){
    Y_len <- rowSums(!is.na(data[,Y_nm])) # find length of functional outcome

    ID <- sapply(1:nrow(data), function(x)
            rep(data[,cluster][x],
                each = Y_len[x]) )
    if(is.list(ID)){
      ID <- do.call(c, ID)
    }else{
      ID <- as.vector(ID)
    }


    if(!is.null(time) & cov.type == "ar1"){
      time.vec <- sapply(1:nrow(data), function(x)
        rep(data[,time][x],
            each = Y_len[x]) )
      time.vec <- as.vector(time.vec)
    }

    colnames(data[[Y_nm]]) <- 1:ncol(data[, Y_nm]) # rename columns -- added 3/20/25 to avoid issues with yindex.vec below
    Y_ <- colnames(data[, Y_nm])

  }else{
    Y_len <- sum(!is.na(data[,Y_nm])) # find length of functional outcome
    ID <- as.vector(data[,cluster])
    time.vec <- as.vector(data[,time])
    Y_ <- Y_nm

    }

  X_ <- colnames(MM)

  if(length(ID) != length(fit_pffr$model[, Y_nm])){
    stop("Some rows in dataset are being dropped in pffr() fit, remove these rows first")
  }

  if(!is.null(time) & cov.type == "ar1"){
    # prepare time for ar1
    dt <- data.table(Y = as.numeric(fit_pffr$model[, Y_nm]),
                       cluster = ID,
                       yindex.vec = fit_pffr$model$yindex.vec,
                       time = time.vec,
                       MM)

    setorder(dt, cluster, time, yindex.vec) # CRITICAL THAT IT IS IN THIS ORDER
    ID <- as.vector(dt$cluster)

    data <- copy(dt)
    rm(dt)
  }
  else{

    data <- data.table(Y = fit_pffr$model[,Y_nm],
                       cluster = ID,
                       yindex.vec = yindex.vec,
                       MM)
  }

  if (any(Y_ %in% X_)) {
    stop(paste("Outcome variable", Y_, "cannot also be a predictor"))
  }

  if(index == "time")   rho.smooth <- FALSE # cannot smooth since there is only one estimate
  namesd <- paste0("d.", X_)
  dx <- copy(data)
  dx[, cname_ := cluster ]

  formula <- stats::update(formula, Y ~ .)

  clusters <- unique(dx[, cname_])
  N_clusters <- length(clusters)

  if(is.null(cv.grid)){
    cv.grid <- vector(length = 3, "list")
    # pre-grid (never include 0 in pre-grid!)
    cv.grid[[1]] <- sort(unique( c(10^(-6:4),
                                   10^(-6:4) * 0.25,
                                   10^(-6:4) * 0.5,
                                   10^(-6:4) * 0.75,
                                   10^(-6:4) * 0.9)))
    cv.grid[[2]] <- c(0, 1e-2, 0.1, 1, 10, 100)
    cv.grid[[3]] <- c(0, 0.5, 0.75, 1, 1.3, 2, 5)

  }

  if(cov.type == "independence"){
    # fast work around to induce working-independence structure
    cov.type <- "ar1"
    weighted <- FALSE
  }else{
    weighted <- TRUE
  }

  if(exact){
    # gaussian exact
    mod.fit <- fun.gee1step.exact(orig.data=data, dx=dx, formula, family, X_, Y_, namesd,
                                 N_clusters, clusters = clusters, yindex.vec, bs=bs,
                                 cov.type = cov.type, glmfit = fit_pffr,
                                 parallel = parallel, sandwich = sandwich, cv = cv,
                                 cv.grid = cv.grid, weighted = weighted,
                                  V.inv = V.inv, time = time, index = index)
  }else{
    # one-step
    mod.fit <- fun.gee1step.dist(orig.data=data, dx=dx, formula, family, X_, Y_, namesd,  bs=bs,
                                 N_clusters, clusters = clusters, yindex.vec, glmfit = fit_pffr,
                                 cov.type = cov.type, parallel = parallel, sandwich = sandwich, cv = cv,
                                 cv.grid = cv.grid, rho.smooth = rho.smooth, weighted = weighted,
                                 return.early = return.early, joint.CI = joint.CI, V.inv = V.inv, time = time, index = index)
  }


  # update model parameters
  mod.fit <- model.update(mod.fit = mod.fit, MM = MM)

  result <- append(mod.fit,
    list(
      call = match.call(),
      formula = orig.formula,
      family = family,
      outcome = Y_,
      xnames = X_,
      model.data = MM,
      data = data,
      pffr_initial.fit = fit_pffr,
      cluster_sizes = as.vector(dx[, .N, keyby = cname_][, N])
    ))
  attr(result, "class") <- "fgee1step"

  return(result)
}



# Sandwich estimator with working independence

#' Estimate parameters using one-step algorithm
#' @importFrom sanic solve_chol
#' @importFrom Rfast rowMaxs ar1 colMaxs
#' @import data.table mgcv refund Matrix SuperGauss MASS
#' @param formula an object of class "formula": a symbolic description of the model to be fitted.
#' @param data a required data frame or data.table containing the variables in the model.
#' @param cluster the name of the field that identifies the clusters.
#' @param family the distribution family: gaussian, binomial, and poisson
#' @param ... currently disregarded
#' @references Lipsitz, S., Fitzmaurice, G., Sinha, D., Hevelone, N., Hu, J.,
#' & Nguyen, L. L. (2017). One-step generalized estimating equations with large
#' cluster sizes. Journal of Computational and Graphical Statistics, 26(3), 734-737.
#' @return a "fgee1step" object
#' @examples
#' geefit <- gee1step(y ~ x1 + x2 + x3, data = sampData_binomial,
#'   cluster = "site", family = "binomial")
#' geefit
#'
#' @export
fun.sandwich.ind <- function(formula, data, cluster, family, beta.hat = NULL,
                         algorithm = "bam", knots = NULL, # smooth_method="fREML",
                         bs = "tp", parallel = FALSE, cov.type = "exchangeable",
                         # wd = "/Users/loewingergc/Desktop/NIMH Research/functional_GEE/fastFGEE/sources/",
                         sandwich = FALSE, cv = TRUE, cv.grid, exact = FALSE, weighted = TRUE, rho.smooth = FALSE,
                          m.pffr = c(2, 1), V.inv = NULL, time = NULL, index = "yindex.vec", ...) {

  # "declare" vars to avoid global NOTE
  # source(paste0(wd, "penalty_matrix_construct.R"))
  # source(paste0(wd, "internals.R"))
  # suppressPackageStartupMessages(Rcpp::sourceCpp(paste0(wd, "srcRcpp.cpp")))

  orig.formula <- formula

  # "declare" vars to avoid global NOTE

  cname_ <- NULL
  Y <- NULL
  N <- NULL

  ### Check arguments

  if ( ! inherits(formula, "formula" ) ) {
    stop("The argument `formula` is not properly specified")
  }

  if  (! is.data.frame(data) ) {
    stop("Data must be a data.frame or data.table")
  }

  if ( ! (cluster %in% names(data)) ) {
    stop(paste("Cluster variable", cluster, "is not in data set"))
  }

  Y_nm <- all.vars(formula)[1]
  L <- ncol(data[, Y_nm])

  if(is.null(knots))  knots <- round(L/4)
  message("Fitting pffr() model")
  # need to order by cluster, yindex.vec, time in the appropriate order
  fit_pffr <- refund::pffr(formula = formula,
                           algorithm = algorithm,
                           family = family,
                           discrete = (algorithm == "bam"),
                           bs.yindex = list(bs = bs,
                                            k = knots,
                                            m = m.pffr),
                           data = data)

  if(!is.null(beta.hat)){
    fit_pffr$coefficients <- beta.hat # use this if feed pre-existing estimates
  }

  yindex.vec <- fit_pffr$model$yindex.vec # functional domain index
  MM <- suppressWarnings(stats::model.matrix(fit_pffr))

  Y <- fit_pffr$model$Y
  Y_len <- rowSums(!is.na(data[,Y_nm])) # find length of functional outcome

  # ID and time
  ID <- sapply(1:nrow(data), function(x)
    rep(data[,cluster][x],
        each = Y_len[x]) )
  ID <- as.vector(ID)
  if(!is.null(time) & cov.type == "ar1"){
    time.vec <- sapply(1:nrow(data), function(x)
      rep(data[,time][x],
          each = Y_len[x]) )
    time.vec <- as.vector(time.vec)
  }

  colnames(data[[Y_nm]]) <- 1:ncol(data[, Y_nm]) # rename columns -- added 3/20/25 to avoid issues with yindex.vec below
  Y_ <- colnames(data[, Y_nm])
  X_ <- colnames(MM)

  if(!is.null(time) & cov.type == "ar1"){

    dt <- data.table(Y = as.numeric(fit_pffr$model[,Y_nm]),
                     cluster = ID,
                     yindex.vec = fit_pffr$model$yindex.vec,
                     time = time.vec,
                     MM)

    setorder(dt, cluster, time, yindex.vec) # CRITICAL THAT IT IS IN THIS ORDER
    ID <- as.vector(dt$cluster)

    data <- copy(dt)
    rm(dt) # time.mat,
  }
  else{

    data <- data.table(Y = fit_pffr$model[,Y_nm],
                       cluster = ID,
                       yindex.vec = yindex.vec,
                       MM)
  }

  if (any(Y_ %in% X_)) {
    stop(paste("Outcome variable", Y_, "cannot also be a predictor"))
  }

  if(index == "time")   rho.smooth <- FALSE # cannot smooth since there is only one estimate
  namesd <- paste0("d.", X_)
  dx <- copy(data)
  dx[, cname_ := cluster ]

  formula <- stats::update(formula, Y ~ .)

  clusters <- unique(dx[, cname_])
  N_clusters <- length(clusters)

  mod.fit <- fun.gee1step.sandwich(orig.data=data, dx=dx, formula, family, X_, Y_, namesd,
                               N_clusters, clusters = clusters, yindex.vec, glmfit = fit_pffr,
                               cov.type = "ar1", parallel = parallel, sandwich = sandwich, cv = cv,
                               cv.grid = cv.grid, rho.smooth = rho.smooth,
                               V.inv = V.inv, time = time, index = index)


  result <- append(mod.fit,
    list(
      call = match.call(),
      formula = orig.formula,
      family = family,
      outcome = Y_,
      xnames = X_,
      model.data = MM,
      data = data,
      pffr_initial.fit = fit_pffr,
      cluster_sizes = as.vector(dx[, .N, keyby = cname_][, N])
    ))
  attr(result, "class") <- "fgee1step"

  return(result)
}

