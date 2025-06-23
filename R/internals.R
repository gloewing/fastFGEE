# Internal functions
D.estimate <- function(dx,
                       namesd,
                       cname_,
                       family,
                       corr = "exchangeable",
                       dt = TRUE,
                       exact = FALSE,
                       index = "yindex.vec" # c("yindex.vec", "time")
                       ){
  # wrapper function for estimate of D depending on correlation structure
  dx <- data.table::copy(dx) # prevent issues with global enviornment

  if(corr == "exchangeable"){
    if(dt){
      # data.table version
      return( .getD_dt(dx, namesd, cname_, family, exact, index = index) )
    }else{
      if(exact){
        return( .get_rhs.exact(dx, namesd, cname_, family, index = index) )
      }else{
        return( .getD(dx, namesd, cname_, family) ) # doesn't work for different index (only yindex.vec)
      }
    }

  }else if(corr == "ar1"){
    return( .getD_AR_dt(dx, namesd, cname_, family, exact, index = index) )
  }
}

# Internal functions
D.estimate.exact <- function(dx,
                               namesd,
                               cname_,
                               family,
                               corr = "exchangeable",
                               dt = TRUE,
                             index = "yindex.vec" #
                             ){
  # wrapper function for estimate of D depending on correlation structure
  dx <- data.table::copy(dx) # prevent issues with global enviornment

  if(corr == "exchangeable"){
    if(dt){
      # data.table version
      return( .getD_exact_dt(dx, namesd, cname_, index = index) )
    }else{
      return( .get_rhs.exact(dx, namesd, cname_) )
    }
  }else if(corr == "ar1"){
    # only data.table version
    return( .getD_AR_exact_dt(dx, namesd, cname_, family, index = index) )
  }
}


W.estimate <- function(dx,
                       namesd,
                       cname_,
                       family,
                       corr = "exchangeable",
                       dt = TRUE,
                       index = "yindex.vec" # c("yindex.vec", "time")
                        ){
  # wrapper function for estimate of D depending on correlation structure
  dx <- data.table::copy(dx) # prevent issues with global enviornment

  if(corr == "exchangeable"){
    if(dt){
      # data.table version
      return( .getW_dt(dx, namesd, cname_, family, index = index) )
    }else{
      return( .getW(dx, namesd, cname_, family) )
    }
  }else if(corr == "ar1"){
    # only data.table version
    return( .getW_AR_dt(dx, namesd, cname_, family, index = index) )
  }
}

# update model object with fitted values etc.
model.update <- function(mod.fit, MM){
  mod.fit$model$coefficients <- mod.fit$beta
  mod.fit$linear.predictors <- MM %*% mod.fit$beta

  # fitted values
  link.fn <- mod.fit$model$family[[2]]
  if(link.fn == "identity"){
    mod.fit$model$fitted.values <- mod.fit$linear.predictors
  }else if(link.fn == "log"){
    mod.fit$model$fitted.values <- exp(mod.fit$linear.predictors)
  }else if(link.fn == "logit"){
    mod.fit$model$fitted.values <- exp(mod.fit$linear.predictors) /
      (1 + exp(mod.fit$linear.predictors))
  }
  return(mod.fit)
}

# joint.basis.np bootstraps in terms of original betas (not f.hat's)
joint.basis.np <- function(glmfit,
                     di,
                     wi,
                     beta,
                     penalty_diag,
                     bs,
                     var.mat,
                     exact = FALSE,
                     grid.size = 500,
                     NN.sim = 5000){

  # calculate qn in terms of basis

  # resample
  boot.samp <- suppressMessages(
                var.est(di = di, wi = wi, beta = beta,
                        penalty_diag = penalty_diag, B = NN.sim,
                        sandwich = "fastBoot", exact = exact,
                        return.boot = TRUE) )$boot

  pp <- length(glmfit$smooth)
  idx.cs <- sapply(1:pp, function(p) ncol(glmfit$smooth[[p]]$D))
  idx.cs <- c(0, cumsum(idx.cs))
  qn <- rep(NA, pp)
  var.vec <- diag(var.mat)
  for(k in 1:pp){
    b.idx <- seq(idx.cs[k]+1, idx.cs[k+1]) # indices for these coefs
    b.mat <- boot.samp[, b.idx] # boot samples for these coefs

    # calculate statistic
    z.beta <- abs( beta[b.idx] - t(b.mat) ) / sqrt(var.vec[b.idx])
    un <- Rfast::colMaxs(z.beta, value = TRUE)
    qn[k] <- stats::quantile(un, 0.95)
  }

  return(qn)
}

# joint.np bootstraps in terms of f.hat's
joint.np <- function(glmfit,
                     di,
                     wi,
                     beta,
                     penalty_diag,
                     bs,
                     var.mat,
                     grid.size = 500,
                     NN.sim = 5000){

  # resample
  boot.samp <- suppressMessages(
                 var.est(di = di, wi = wi, beta = beta,
                        penalty_diag = penalty_diag, B = NN.sim,
                        sandwich = "fastBoot", exact = FALSE,
                        return.boot = TRUE) )
  var.mat <- boot.samp$cov
  boot.samp <- boot.samp$boot

  # grid size is at least as big as original grid of points
  yind <- glmfit$pffr$yind
  if(length(yind) > grid.size){
    grid.size <- length(yind)
  }
  x <- seq(min(yind), max(yind), length = grid.size) # grid of points for functional domain -- can alter this if need different grid
  pp <- length(glmfit$smooth)
  idx <- sapply(1:pp, function(p) ncol(glmfit$smooth[[p]]$D))
  idx.cs <- c(0,cumsum(idx))
  grid.idx <- c(0, cumsum(rep(grid.size, pp)))
  qn <- rep(NA, pp)
  sm.X.ls <- lapply(1:pp, function(kk)
                    mgcv::smoothCon( s(x, bs=bs, k = idx[kk]), data=data.frame(x) )[[1]]$X ) # design matrix
  XX <- Matrix::bdiag(sm.X.ls)
  var.vec <- rowSums( XX * t(var.mat %*% t(XX)) )
  rm(XX)
  for(k in 1:pp){
    b.idx <- seq(idx.cs[k]+1, idx.cs[k+1]) # indices for these coefs
    b.mat <- boot.samp[, b.idx] # boot samples for these coefs
    sm.X <- sm.X.ls[[k]]
    f.boot.hat <- tcrossprod(sm.X, b.mat) # calculate f.hat for bootstraped betas

    # real f.hat
    grid.idx.k <- seq(grid.idx[k]+1, grid.idx[k+1]) # indices for f.hat
    f.hat <- as.numeric( sm.X %*% beta[b.idx] )

    # calculate statistic
    z.beta <- abs(f.hat - f.boot.hat) / sqrt(var.vec[grid.idx.k])
    un <- Rfast::colMaxs(z.beta, value = TRUE)
    qn[k] <- stats::quantile(un, 0.95)
  }

  return(qn)
}


joint.qn <- function(glmfit, NN.sim = 100000){

  p.vec <- sapply(glmfit$smooth, function(xx) ncol(xx$S[[1]])) # dimensions of smooths
  p.vec[1] <- p.vec[1] + 1 # intercept smooth has extra dimension for non-penalized intercept
  p.vec <- c(0, cumsum(p.vec))

  qn <- rep(0, length = length(glmfit$smooth))

  for(i in 1:length(qn)) {
    cov.idx <- (p.vec[i]+1):(p.vec[i+1])
    Sigma <- glmfit$Vp[cov.idx, cov.idx]
    sqrt_Sigma <- sqrt(diag(Sigma))
    S_scl <- Matrix::Diagonal(x = 1 / sqrt_Sigma)
    Sigma <- as.matrix(S_scl %*% Sigma %*% S_scl)
    zero_vec <- rep(0, nrow(Sigma))
    x_sample <- abs(MASS::mvrnorm(NN.sim, zero_vec, Sigma))
    un <- Rfast::rowMaxs(x_sample, value = TRUE)
    qn[i] <- stats::quantile(un, 0.95)
  }

  return(qn)
}

corr.estimate <- function(dx,
                          corr = "exchangeable",
                          rho.smooth = FALSE,
                          ar = "yw", # c("mom", "yw")
                          glmfit = NULL,
                          index = "yindex.vec"){

  dx <- data.table::copy(dx) # prevent issues with global enviornment
  dx[, index.vec := get(index)]

  if(corr == "exchangeable"){
    ### Estimate ICC
    drho <- dx[,
               list(
                 .N,
                 sum_r = sum( resid ),
                 uss_r = sum( resid^2 )
               ), keyby = .(cname_, index.vec)]

    drho[, wt_ij := ( N * (N-1) / 2), by = index.vec]
    drho[, rho_ij := sum_r^2 - uss_r, by = index.vec ]
    rho <- drho[, (sum(rho_ij)/2) / sum(wt_ij), by = index.vec]

    if(rho.smooth){
      drho[, wt_ij := (N*(N-1) / 2), keyby = .(cname_, index.vec)]
      drho[, rho_ij := sum_r^2 - uss_r, keyby = .(cname_, index.vec) ]
      rho.knots <- glmfit$smooth[[1]]$bs.dim + 1
      drho$rho_fits <- mgcv::gam(rho_ij ~ s(index.vec, k = rho.knots),
                                 data = drho)$fitted.values
      rho <- drho[, (sum(rho_fits)/2) / sum(wt_ij), keyby = .(index.vec)]
    }

    rho[, V1 := ifelse(abs(V1) > 0.999, 0.999 * sign(V1), V1)]

  }else if(corr == "ar1"){
    ### Estimate ICC
    n.vec <- dx[,.N, keyby = .(cname_, index.vec)]
    if(ar == "yw"){
      drho <- dx[, Rfast::ar1(resid, method = "yw")[2], # method = "yw"
                 keyby = .(cname_, index.vec)][n.vec,
                 on = .(cname_, index.vec)][,
                 n.clust := length(unique(cname_))]

    }else if(ar == "mom"){
      drho <- dx[, as.numeric(resid[-1]) %*% as.numeric(resid[-.N]),
                 keyby = .(cname_, index.vec)][n.vec,
                 on = .(cname_, index.vec)][,
                 n.clust := length(unique(cname_))]
    }

    if(anyNA(drho$V1)){
      message(paste("% missing in corr. estimation:",
                    sum(is.na(drho$V1))/nrow(drho)))
      # missing values
      drho <- na.omit(drho, cols="V1")
      drho[, n.clust := .N, keyby = .(index.vec)]
      drho[, N := .N, keyby = .(cname_)]
    }

    if(rho.smooth){
      rho.knots <- glmfit$smooth[[1]]$bs.dim + 1
      drho$V1 <- mgcv::gam(V1 ~ s(index.vec, k = rho.knots),
                                 data = drho)$fitted.values
    }

    rho <- drho[, sum(V1), by = index.vec] #
    drho <- drho[,c("cname_", "index.vec", "N", "n.clust")][,
                  .SD[1], by = index.vec]
    rho <- rho[drho[,c("index.vec", "N", "n.clust")], on = "index.vec"]

    if(ar == "yw"){
      rho[, V1 := V1 / n.clust]
    }else if(ar == "mom"){
      rho[, V1 := V1 / (n.clust*(N-1))]
    }

    if(index == "time"){
      # sum across time and make all equal
      rho[,V1 := sum(V1) / max(index.vec)]
    }

    rho[, V1 := ifelse(abs(V1) > 0.999, 0.999 * sign(V1), V1)]
  }

  colnames(rho)[2] <- "rho"
  return(rho)
}

var.est <- function(di, wi, beta,
                    penalty_diag, B = 500,
                    sandwich = TRUE, exact = FALSE,
                    clust.vec = NULL,
                    return.boot = FALSE){

  N_clusters <- length(wi)
  clusters <- 1:N_clusters
  if(!is.null(clust.vec)){
    # adjust for varying cluster sizes
    clust.ls <- sapply(clusters, function(ii) length(clust.vec == ii))
    NN <- sum(clust.ls)
  }else{
    clust.ls <- rep(1, length = length(clusters))
    NN <- sum(clust.ls)
  }

  if(exact){
    if(sandwich %in% c(TRUE, FALSE) ){
      di <- lapply(di, function(i) tcrossprod(i) )
      di <- Reduce("+", di) #/N_clusters
      wi <- sanic::solve_chol( Reduce("+", wi) + penalty_diag )
      return( wi %*% di %*% wi )

    }else if(sandwich == "boot"){
      message("bootstrap variance")
      boot.mat <- matrix(NA, ncol = length(beta), nrow = B)
      for(b in 1:B){
        message(paste("bootstrap var iter:",b))
        # draw sample of clusters
        samp <- sample(clusters, N_clusters, replace = TRUE)
        nk <- sum(clust.ls[samp]) / NN # rows of held out set
        W <- Reduce("+", wi[samp])
        D <- Reduce("+", di[samp])
        boot.mat[b,] <- sanic::solve_chol(a = W + nk*penalty_diag,
                                          b = D)  # coefficient estimate
      }

      if(return.boot){
        return(list(cov = cov(boot.mat), boot = boot.mat) )
      }else{
        return( cov(boot.mat) )
      }
    }else if(sandwich == "fastBoot"){
      message("Fast Bootstrap Variance")
      nn <- N_clusters
      boot.mat <- matrix(NA, ncol = length(beta), nrow = B)
      wwi <- sanic::solve_chol( Reduce("+", wi) + penalty_diag )

      for(b in 1:B){
        samp <- sample(clusters, nn, replace = TRUE) # draw sample of clusters
        nk <- NN / sum(clust.ls[samp]) # rows of held out set
        boot.mat[b,] <- wwi %*% Reduce("+", di[samp]) * nk # coefficient estimate
      }

      if(return.boot){
        return(list(cov = cov(boot.mat), boot = boot.mat) )
      }else{
        return( cov(boot.mat) )
      }

      }
    # end exact -----------------------------

    }else{
      # one-step
      if(sandwich == TRUE){
        message("Sandwich Variance Estimator")
        penalty_vec <- as.numeric(penalty_diag %*% beta * N_clusters)
        W <- Reduce("+", wi)/N_clusters
        D <- do.call(cbind, di)
        inf_fn <- sanic::solve_chol(a = W + penalty_diag,
                                    b = D - penalty_vec)
        inf_fn <- inf_fn - rowMeans(inf_fn)
        return( 1/(N_clusters*(N_clusters-1)) * tcrossprod(inf_fn) )
      }else if(sandwich == "boot"){
        message("Bootstrap Variance")
        penalty_vec <- as.numeric(penalty_diag %*% beta * N_clusters)
        boot.mat <- matrix(NA, ncol = length(beta), nrow = B)
        for(b in 1:B){
          # draw sample of clusters
          samp <- sample(clusters, N_clusters, replace = TRUE)
          nk <- sum(clust.ls[samp]) / NN # rows of held out set
          # # -----------------------------------------------
          W <- Reduce("+", wi[samp])/N_clusters
          D <- Reduce("+", di[samp])
          inf_fn <- sanic::solve_chol(a = W + penalty_diag*nk,
                                      b = D - penalty_vec*nk)  # coefficient estimate
          boot.mat[b,] <- beta + inf_fn / N_clusters # coefficient estimate
        }
        return( stats::cov(boot.mat) )
      }else if(sandwich == "fastBoot"){
        message("Fast Bootstrap Variance")
        boot.mat <- matrix(NA, ncol = length(beta), nrow = B)
        wwi <- sanic::solve_chol( Reduce("+", wi)/N_clusters + penalty_diag )
        penalty_vec <- as.numeric(penalty_diag %*% beta * N_clusters)
        for(b in 1:B){
          samp <- sample(clusters, N_clusters, replace = TRUE) # draw sample of clusters
          nk <- NN / sum(clust.ls[samp]) # rows of held out set
          inf_fn <- wwi %*% (Reduce("+", di[samp]) * nk - penalty_vec)
          boot.mat[b,] <- beta + inf_fn / N_clusters # coefficient estimate
        }
        if(return.boot){
          return(list(cov = cov(boot.mat), boot = boot.mat) )
        }else{
          return( cov(boot.mat) )
        }

      }
  }
}


.getW_ar1 <- function(dd, namesd, family = "gaussian") {
  v <- NULL
  rho <- NULL
  time <- NULL
  sqrtv <- NULL

  dd[, sqrtv := sqrt(v)]
  pp <- length(namesd)
  rho <- unique(dd$rho)
  d <- Matrix::Diagonal(x = 1/as.numeric(dd$sqrtv)) %*% as.matrix(dd[, namesd, with = FALSE])

  return( crossprod(d,
                    sanic::solve_chol(a=toeplitz(rho^(dd$time-1)), b=d)) ) # GL added in /(1-rho). Also changed from:  d1 - xrho * d2 %*% t(d2)
}

.getD_ar1 <- function(dd, namesd, family = "gaussian") {
  v <- NULL
  time <- NULL
  rho <- NULL
  sqrtv <- NULL
  dd <- data.table::copy(dd) # prevent issues with global enviornment
  #setorder(dd, time)
  dd[, sqrtv := sqrt(v)]
  dd[, YY := (Y-p) / sqrtv]
  pp <- length(namesd)
  rho <- unique(dd$rho)
  d <- Matrix::Diagonal(x = 1/as.numeric(dd$sqrtv)) %*% as.matrix(dd[, namesd, with = FALSE])

  return( crossprod(d,
                    sanic::solve_chol(a=toeplitz(rho^(dd$time-1)), b=as.numeric(dd$YY)) )) # GL added in /(1-rho). Also changed from:  d1 - xrho * d2 %*% t(d2)
}

.getW_dt <- function(dd, namesd, cname_, family = "gaussian",
                     index = "yindex.vec", second.deriv = FALSE) {
  dd <- data.table::copy(dd) # prevent issues with global enviornment
  dd[, index.vec := get(index)]
  v <- NULL
  pp <- length(namesd)
  dd[, ':=' (resid = (Y - p), sqrtv = sqrt(v))][,
       nn := .N, by = .(cname_, index.vec)]

  # This does the scaling of the design matrix for non-gaussians
  if(family == "poisson"){
    dd3 <- copy(dd)
    dd <- dd[, scl2 := p][, # first derivative of g(E[Y | X]) ignoring design matrix
               sqrtv := sqrtv / scl2] # scaling sqrtv is trick to allow use of code below
  }else if(family == "binomial"){
    dd3 <- copy(dd)
    dd <- dd[, scl2 := p*(1-p)][, # first derivative of g^-1(E[Y | X]) ignoring design matrix
               sqrtv := sqrtv / scl2] # scaling sqrtv is trick to allow use of code below
  }

  d2 <- copy(dd) # copy after adjustments to sqrtv above
  res <- dd[, xrho := rho / ( (1-rho) + nn * rho ),
            by = .(cname_, index.vec)][,.SD[1],
            .SDcols = c("xrho", "rho"),
            by = .(cname_, index.vec)] # removes all the repeats

  d2 <- d2[, (namesd) := lapply(.SD, function(xx) sum(xx/ sqrtv)),
             .SDcols = namesd,
             by = .(cname_, index.vec)][,.SD[1],
              .SDcols = namesd,
              by = .(cname_, index.vec)]

    # tcrossprod(d2) -- each row gives same as tcrossprod(dv(d, sqrtv)) # WORKS FOR BINARY
  d2 <- d2[, crossprod(as.matrix(.SD)),
             .SDcols = namesd,
             by = .(cname_, index.vec)]
    #--------------------------------------------
    # divide all columns by sqrtv
    dd <- dd[, (namesd) := lapply(.SD, function(xx) xx / sqrtv),
             .SDcols = namesd]

    # this should be same as crossprod(X/sqrtv) = ddv(d, v)[1:5,1:5]
    dd <- dd[, crossprod(as.matrix(.SD)),
             .SDcols = namesd,
             by = .(cname_, index.vec)][,
              # add matrix indices to join with d2 below
             # m_idx := rep(1:(pp^2), .N/(pp^2)), # added in to join below
             m_idx :=1:(pp^2), # added in to join below
             by = .(cname_, index.vec)][,
             D1 := V1][,
             V1 := NULL]
    d2 <- d2[res, on = .(cname_, index.vec)][,
              D2 := V1 * xrho][, # scale the d2 values by xrho and delete old column
              V1 := NULL][,
              # m_idx := rep(1:(pp^2), .N/(pp^2)), # add in matrix indices to join with dd
              m_idx := 1:(pp^2), # add in matrix indices to join with dd (should be same as above)
                        by = .(cname_, index.vec)][dd,
              on = .(cname_, index.vec, m_idx)][, # join dd and d2
              ':=' (m = (D1 - D2) / (1 - rho), # this would be negative of Gaussian case so subtract below
                    m_idx = rep(1:(pp^2), .N/(pp^2)) ) ][,
              sum(m), keyby = .(cname_, m_idx)] # sum over yindex

    rm(dd, res)
    clusts <- unique(d2$cname_)
    # ------------------------------------------
    # -------------------------
    # second derivative part
    # -------------------------
    if( (family == "poisson" | family == "binomial") & second.deriv){
      # scl1 <- as.numeric(m_res * (1 - m_res) * (1 - 2*m_res) * Vinv %*% resid)
      # term1 <- crossprod(d * scl1, d)

      # second derivative of g^{-1}(E[Y | X])
      if(family == "binomial"){
        dd3 <- dd3[,glm.scl :=  p * (1 - p) * (1 - 2*p)]
      }else if(family == "poisson"){
        dd3 <- dd3[,glm.scl :=  p]
        }

      dd3 <- dd3[, ':=' (xrho = rho / ((1-rho)*( (1-rho) + nn * rho )),
                         resid = resid)][,
                         r.sum := sum(resid), by = .(cname_, yindex.vec)][,
      # this should be: scl1 <- as.numeric(m_res * (1 - m_res) * (1 - 2*m_res) * Vinv %*% resid)
                    ':=' (T1 = 1/((1-rho)*v) * resid,
                         T2 = xrho/v * r.sum)][,
                         term1_scl1 := glm.scl * (T1 - T2)][,.SD,
                         .SDcols=c("cname_", "yindex.vec", "term1_scl1", namesd)][,
      # should be same as crossprod(X * term1_scl1, X)
                  crossprod(as.matrix(.SD)*term1_scl1, as.matrix(.SD)),
                  .SDcols = namesd,
                  by = .(cname_, yindex.vec)][,
                  m_idx := rep(1:(pp^2), .N/(pp^2))][,
                  sum(V1), keyby = .(cname_, m_idx)][,# sum over yindex.vec
                  term1 := V1 ][,
                  V1 := NULL][d2,  # rename term1 to avoid confusion on join below
                # join term1 and 2 so can calculate difference
                  on = .(cname_, m_idx)][,
                  res := -term1 + V1][,.SD,
                  .SDcols=c("cname_", "res")]
      message("2nd derivative")

      rm(d2)
      #---------------------------------------------------
      return( lapply(clusts, function(ii)
        matrix(dd3[cname_ == ii, res],
               byrow = TRUE,
               nrow = pp, ncol = pp) ) )
    }else if(family == "gaussian" | !second.deriv ){
      return(lapply(clusts, function(ii)
        matrix(d2[cname_ == ii, V1],
               byrow = TRUE,
               nrow = pp, ncol = pp) ) )
    }

}

#############################

.getW <- function(dd, namesd, rho, family = "gaussian") {
  # "declare" vars to avoid global NOTE
  v <- NULL

  ###
  d <- as.matrix(dd[, namesd, with = FALSE])

  v <- dd[, v]
  sqrtv <- dd[, sqrt(v)]
  n <- length(v)
  # if(family != "gaussian"){
  #   # use this as a work around for binomial/Poisson (may be different for
  #   # GLM families if a(\phi) != 1) to make Hessian = t(X) %*% W %*% X,
  #   # where W is diagonal with (see GLM notes) log(p) (for Poisson), or p*(1-p) (for binomial)
  #   d <- d/v
  #   #browser()
  # }

  if(family == "binomial"){
    m_res <- dd[, p]
    Y <- as.numeric(as.matrix(dd[,"Y"]))

    resid <- as.numeric(Y - m_res)
    Vinv <- ((1/(1-rho)) * diag(n) - rho/((1-rho) *
                    (1-rho + n*rho)) * matrix(1, n, n)) / v
    # term 1 (2nd deriv part)
    scl1 <- as.numeric(m_res * (1 - m_res) * (1 - 2*m_res) * Vinv %*% resid)
    term1 <- crossprod(d * scl1, d)

    # term 2  -D^T W D
    scl2 <- m_res * (1 - m_res)  # 1st deriv
    term2 <- -crossprod(d * scl2, Vinv %*% d * scl2)
    return( term1 + term2 )
    # ---------------------------------------
  }else if(family == "gaussian"){
    xrho <- rho / ( (1-rho) + nrow(d) * rho )
    d1 <- ddv(d, v)
    d2 <- dv(d, sqrtv)

    # return( (d1 - xrho * tcrossprod(d2)) )
    return( (d1 - xrho * tcrossprod(d2)) / (1-rho) )

  }

}

##########################
.getW_true <- function(dd, namesd, V.inv, rho = NULL, family = "gaussian") {
    d <- as.matrix(dd[order(dd$yindex.vec), namesd, with = FALSE])
    return( crossprod(d, V.inv %*% d) ) # GL added in /(1-rho). Also changed from:  d1 - xrho * d2 %*% t(d2)
}

.getD_true <- function(dd, namesd, V.inv, rho = NULL, family = "gaussian") {
  # only works for GLS/exact solution NOT one-step (one-step would use residual instead of Y)
  XX <- as.matrix(dd[order(dd$yindex.vec), namesd, with = FALSE])

  ## commented out Oct 1, 2024
  dd[, r := (Y - p) ] #
  rr <- as.numeric(dd$r)

  return( crossprod(XX, V.inv %*% rr) )
}

.getD_true_exact <- function(dd, namesd, V.inv, rho = NULL, family = "gaussian") {
  # only works for GLS/exact solution NOT one-step (one-step would use residual instead of Y)
  XX <- as.matrix(dd[order(dd$yindex.vec), namesd, with = FALSE])

  # dd[, r := (Y - p) / v ] #
  YY <- as.numeric(dd$Y)

  return( crossprod(XX, V.inv %*% YY) )
}


.get_rhs <- .getD <- function(dd, namesd, rho, family = "gaussian") {
  X <- as.matrix(dd[,namesd, with = FALSE])
  # added below in 12/20/24 to account for removing from fun.gee1step.dist
  if(family == "binomial"){
    X <- X * dd[, p*(1-p)]
  }else if(family == "poisson"){
    X <- X * dd[, p]
  }
  ## commented out Oct 1, 2024
  dd[, r := (Y - p) / v ] #
  r <- as.numeric(dd$r)
  n <- nrow(X)

  scl1 <- 1/(1-rho)
  scl2 <- rho/((1-rho) * (1-rho + n*rho))

  return( scl1 * crossprod(X, r) - scl2 * colSums(X) * sum(r) )
  #return( crossprod(X, V %*% r)  ) # equivalent for binomial case and guassian 12-20-24
}


.getW_AR_dt <- function(dd, namesd, cname_, family = "gaussian",
                        index = "yindex.vec", second.deriv = FALSE) {

  dd <- data.table::copy(dd) # prevent issues with global enviornment
  dd[, index.vec := get(index)]

  if(index == "time"){
    dd[, ar.idx := yindex.vec] # ensures AR1 structure across functional domain
  }else if(index == "yindex.vec"){
    dd[, ar.idx := time] # ensures AR1 structure across timepoints (longitudinal)
  }

  v <- NULL
  pp <- length(namesd)
  clusts <- unique(dd$cname_)
  dd[, sqrtv := sqrt(v)][,
       nn := .N, by = .(cname_, index.vec)]
  # [,
  #     # divide all columns of design mat by sqrtv so that both A^{-1/2} terms in V^{-1} are incorporated
  #     # important for both term1 and term2
  #      (namesd) := lapply(.SD, function(xx) xx / sqrtv),
  #      .SDcols = namesd]

  # scl2 does the scaling of the design matrix for non-gaussians to take scaling
  # from D = d/(d \beta) \mu into design matrix X
  # dividing by sqrtv so that both A^{-1/2} terms in V^{-1} are incorporated
  # important for both term1 and term2
  if(family == "poisson"){
    dd2 <- copy(dd)
    dd <- dd[, scl2 := p / sqrtv]
  }else if(family == "binomial"){
    dd2 <- copy(dd)
    dd <- dd[, scl2 := p*(1-p) / sqrtv]
  }else if(family == "gaussian"){
    dd <- dd[, scl2 := 1 / sqrtv]
  }

  # added in 3-23-25 from below
  # divide all columns of design mat by sqrtv so that both A^{-1/2} terms in V^{-1} are incorporated
  # important for both term1 and term2
  dd[, (namesd) := lapply(.SD, function(xx) xx * scl2),
     .SDcols = namesd]

  # ------------------------------------------
  # term 1  -D^T W D
  #scl1 <- m_res * (1 - m_res)  # 1st deriv
  #term1 <- -crossprod(d * scl1, Vinv %*% d * scl1)
  # ------------------------------------------
  #--------------------------------------------
  ## AR1 components
  #ar.exp := seq(0, .N-1), # calculate exponent to raise each rho to
  # by = .(cname_, index.vec)][,
  term1 <- dd[, ar.acf := rho^(ar.idx-1)][, # calculate rho
              crossprod(
                  as.matrix(.SD), # * scl2   -- (swapped out3-23-25) scl2 deals scales X for non-gaussian D term ( d/(d/beta) \mu )
                  SuperGauss::Toeplitz$new(N = .N,
                                           acf = ar.acf)$solve(as.matrix(.SD),
                                                               method = "gschur", tol = 1e-8) # solve system of linear equations
                        ),
                .SDcols = namesd, by = .(cname_, index.vec)][,
                m_idx := 1:(pp^2), by = .(cname_, index.vec)][,
                sum(V1), by = .(cname_, m_idx)][, # sum over yindex
                T1 := -V1][, # rename and make negative (so we add T1+T2 below)
                V1 := NULL]

  rm(dd)
  # end term 1
  # ------------------------------------------

  # ------------------------------------------
  # term 2 (2nd deriv part)
  # Vinv <- ((1/(1-rho)) * diag(n) - rho/((1-rho) *
  #           (1-rho + n*rho)) * matrix(1, n, n)) / v
  # # term 2 (2nd deriv part)
  # scl2 <- m_res * (1 - m_res) * (1 - 2*m_res) * Vinv %*% resid
  # term2 <- crossprod(d * as.numeric(scl2),  d)
  # ------------------------------------------

  if( (family == "poisson" | family == "binomial") & second.deriv){
    # scl1 <- as.numeric(m_res * (1 - m_res) * (1 - 2*m_res) * Vinv %*% resid)
    # term1 <- crossprod(d * scl1, d)

    # second derivative of g^{-1}(E[Y | X])
    if(family == "binomial"){
      dd2 <- dd2[,glm.scl :=  p * (1 - p) * (1 - 2*p)]
    }else if(family == "poisson"){
      dd2 <- dd2[,glm.scl :=  p]
    }

    term2 <- dd2[, residuals := (Y - p)][,
                ar.exp := ar.idx-1, # calculate exponent to raise each rho to
                by = .(cname_, index.vec)][,
                ar.acf := rho^ar.exp][, # calculate rho
                Vinv.r := SuperGauss::Toeplitz$new(N = .N, # solve system of linear equations
                                                   acf = ar.acf)$solve(residuals,
                                                                       method = "gschur", tol = 1e-8),
                by = .(cname_, index.vec)][,
                Vinv.r := glm.scl * Vinv.r][, # scale by second derivative terms glm.scl
                crossprod(as.matrix(.SD),
                          as.matrix(.SD) * Vinv.r), # calculate t(X) %*% V.inv %*% r
                .SDcols = namesd, by = .(cname_, index.vec)][,
                m_idx := 1:(pp^2), by = .(cname_, index.vec)][, # matrix indices
                sum(V1), by = .(cname_, m_idx)][,# sum over index.vec
                T2 := V1 ][, # rename to avoid confusion
                V1 := NULL][term1, # join term1 and 2 so can calculate difference
                on = .(cname_, m_idx)][,
                # res := T2 + T1][,.SD,  # add them here because we made T1 negative above
                 res := -(T2 + T1)][,.SD,  # seems to work by CIs appear small: reversed sign 3-21-25
                #res := (T2 - T1)][,.SD,  # cannot tell if this is correct or above
                .SDcols=c("cname_", "res")]
    rm(dd2)
    #---------------------------------------------------
    return(
      # lapply(clusts, function(ii)
      #   matrix(-term1[cname_ == ii, T1],  # negative is correct: this negative was added on 3-3-25!!!
      #          byrow = TRUE,
      #          nrow = pp, ncol = pp) ) )
      lapply(clusts, function(ii)
      matrix(term2[cname_ == ii, res],
             byrow = TRUE,
             nrow = pp, ncol = pp) )
    )
  }else if(family == "gaussian" | !second.deriv ){
    return(lapply(clusts, function(ii)
      matrix(-term1[cname_ == ii, T1],  # negative is correct: this negative was added on 3-3-25!!!
             byrow = TRUE,
             nrow = pp, ncol = pp) ) )
  }

}

.getD_AR_dt <- function(dd, namesd, cname_, family = "gaussian",
                        exact = FALSE, index = "yindex.vec") {

  d <- data.table::copy(dd) # prevent issues with global enviornment
  d[, index.vec := get(index)]

  if(index == "time"){
    d[, ar.idx := yindex.vec] # ensures AR1 structure across functional domain
  }else if(index == "yindex.vec"){
    d[, ar.idx := time] # ensures AR1 structure across timepoints (longitudinal)
  }

  res <- NULL
  residuals <- NULL
  clusts <- unique(d$cname_)
  pp <- length(namesd)

  # added below in 12/20/24 to account for removing from fun.gee1step.dist
  if(family == "binomial"){
    d[,namesd] <- d[,namesd, with = FALSE] * d[, p*(1-p)]
  }else if(family == "poisson"){
    d[,namesd] <- d[,namesd, with = FALSE] * d[, p]
  }

  ## AR1 components
  if(exact){
    res <- d[, sqrtv := sqrt(v)][,
               residuals := Y / sqrtv] # No residuals in exact solution
  }else{
    res <- d[, sqrtv := sqrt(v)][,
                residuals := (Y - p) / sqrtv ]
  }
  # res[, ar.exp := seq(0, .N-1), # calculate exponent to raise each rho to
  #by = .(cname_, yindex.vec)]
  res <- res[, ar.acf := rho^(ar.idx-1)][, # calculate rho
             Vinv.r := SuperGauss::Toeplitz$new(N = .N, # solve system of linear equations
                                                acf = ar.acf)$solve(residuals,
                                                                    method = "gschur", tol = 1e-8),
             by = .(cname_, index.vec)][,
             Vinv.r := Vinv.r / sqrtv][, # scale by sd
                crossprod(as.matrix(.SD), Vinv.r), # calculate t(X) %*% V.inv %*% r
                .SDcols = namesd, by = .(cname_, index.vec)][,
                m_idx := 1:pp, by = .(cname_, index.vec)][, # rep(1:(pp^2), .N/(pp^2))
                sum(V1), by = .(cname_, m_idx)] # sum over yindex

  # rm(dd)
  #
  return(lapply(clusts, function(ii) as.numeric(res[cname_ == ii, V1]) ) )
}

.getD_AR_exact_dt <- function(dd, namesd, cname_,
                              family = "gaussian", index = "yindex.vec") {
  # identical to non-exact but family uncessary and residuals are Y/sqrtv instead
  # of (Y-p)/sqrtv

  dd <- data.table::copy(dd) # prevent issues with global enviornment
  dd[, index.vec := get(index)]

  if(index == "time"){
    dd[, ar.idx := yindex.vec] # ensures AR1 structure across functional domain
  }else if(index == "yindex.vec"){
    dd[, ar.idx := time] # ensures AR1 structure across timepoints (longitudinal)
  }

  clusts <- unique(dd$cname_)
  pp <- length(namesd)
  ## AR1 components
  res <- dd[, sqrtv := sqrt(v)][,
              residuals := Y / sqrtv ][, #exact version uses Y/sqrtv instead of (Y-p)/sqrtv
              ar.exp := ar.idx-1, # calculate exponent to raise each rho to
              by = .(cname_, index.vec)][,
              ar.acf := rho^ar.exp][, # calculate rho
              Vinv.r := SuperGauss::Toeplitz$new(N = .N, # solve system of linear equations
                                                 acf = ar.acf)$solve(residuals,
                                                                     method = "gschur", tol = 1e-8),
              by = .(cname_, index.vec)][,
              Vinv.r := Vinv.r / sqrtv][, # scale by sd
              crossprod(as.matrix(.SD), Vinv.r), # calculate t(X) %*% V.inv %*% r
              .SDcols = namesd, by = .(cname_, index.vec)][,
              m_idx := 1:pp, by = .(cname_, index.vec)][, # rep(1:(pp^2), .N/(pp^2))
              sum(V1), by = .(cname_, m_idx)] # sum over yindex

  rm(dd)

  return(lapply(clusts, function(ii) as.numeric(res[cname_ == ii, V1]) ) )
}

# exangeable D for exact
.getD_exact_dt <- function(dd, namesd, cname_, index = "yindex.vec") {

  dd <- data.table::copy(dd) # prevent issues with global enviornment
  dd[, index.vec := get(index)]

  ## commented out Oct 1, 2024
  dd[, residuals := Y / v ] # exact version is same but Y/v instead of (Y-p)/v
  dd[, ':=' (resid.sum = sum(residuals),
             nn = .N),
     by = .(cname_, index.vec)]

  dd[, ':=' (scl1 = 1/(1-rho),
             scl2 = rho/((1-rho) * (1-rho + nn*rho)) ),
     by = .(cname_, index.vec)]

  # 3. Multiple columns in place
  dd <- dd[, (namesd) := lapply(.SD, function(xx) sum(scl1*residuals * xx) -
                                  sum(scl2 * resid.sum * xx)),
           .SDcols = namesd,
           by = .(cname_, index.vec)]

  dd <- dd[,.SD[1], .SDcols = namesd, by = .(cname_, index.vec)] # removes all the repeats

  dd <- dd[, lapply(.SD, sum), .SDcols = namesd, by = cname_]

  return(lapply(1:nrow(dd), function(xx) as.numeric(dd[xx, ..namesd]) ))
}

.getD_dt <- function(dd, namesd, cname_, family = "gaussian",
                     exact = FALSE, index = "yindex.vec") {

  dd <- data.table::copy(dd) # prevent issues with global enviornment
  dd[, index.vec := get(index)]

  # added below in 12/20/24 to account for removing from fun.gee1step.dist
  if(family == "binomial"){
    dd[,namesd] <- dd[,namesd, with = FALSE] * dd[, p*(1-p)]
  }else if(family == "poisson"){
    dd[,namesd] <- dd[,namesd, with = FALSE] * dd[, p]
  }

  ## commented out Oct 1, 2024
  if(exact){
    dd[, residuals := Y / v ] # exact version is same but Y/v instead of (Y-p)/v
  }else{
    dd[, residuals := (Y - p) / v ]
  }
  dd[, ':=' (resid.sum = sum(residuals),
             nn = .N),
     by = .(cname_, index.vec)]

  dd[, ':=' (scl1 = 1/(1-rho),
             scl2 = rho/((1-rho) * (1-rho + nn*rho)) ),
     by = .(cname_, index.vec)]

  # 3. Multiple columns in place
  dd <- dd[, (namesd) := lapply(.SD, function(xx) sum(scl1*residuals * xx) -
                                  sum(scl2 * resid.sum * xx)),
           .SDcols = namesd,
           by = .(cname_, index.vec)]

  dd <- dd[,.SD[1], .SDcols = namesd, by = .(cname_, index.vec)] # removes all the repeats

  dd <- dd[, lapply(.SD, sum), .SDcols = namesd, by = cname_]

  return(lapply(1:nrow(dd), function(xx) as.numeric(dd[xx, ..namesd]) ))
}

.get_rhs.exact <- function(dd, namesd, rho, family = "gaussian") {

  X <- as.matrix(dd[,namesd, with = FALSE])
  ## commented out Oct 1, 2024
  dd[, r := Y / v ] #
  r <- as.numeric(dd$r)
  n <- nrow(X)

  scl1 <- 1/(1-rho)
  scl2 <- rho/((1-rho) * (1-rho + n*rho))
  return( scl1 * crossprod(X, r) - scl2 * colSums(X) * sum(r) )

  # they are the same for exact: confirmed 12-20-24
  # cbind(scl1 * crossprod(X, r) - scl2 * colSums(X) * sum(r),
  #       crossprod(X, solve(V) %*% Y))
}

# Internal functions

.getLHS <- function(dd, namesd) {

  # "declare" vars to avoid global NOTE
  v <- NULL
  resid <- NULL

  ###
  d <- as.matrix(dd[, namesd, with = FALSE])
  sqrtv <- dd[, sqrt(v)]
  resid <- dd[, resid]

  return( sum( crossprod(d, resid / sqrtv) ) ) # resid is already divided by sqrt(v) so this should be equivalent to: (Y - \mu) / v

}

# construct penalty matrix
penalty.construct <- function(model, lambda = NULL, nonSmooths = 1){

  # penalty
  S <- model$smooth
  L <- length(S)
  if(is.null(lambda))  lambda <- model$sp
  lambda <- as.numeric(lambda)

  dum_pen <- matrix(0, ncol = nonSmooths+1, nrow = nonSmooths+1) # dummy terms
  # extract penalty matrices and scale by lambda -- verified (at least unweighted case)
  pen_mat <- lapply(1:L, function(i) S[[i]]$S[[1]] * lambda[i] ) # * (model$sp[i]) / (S[[i]]$S.scale) #  * (model$sp[i]) / S[[i]]$S.scale # * sqrt(model$sp[i])
  pen_mat <- Matrix::bdiag(pen_mat) # concatenate matrices into block diagonal

  return( as.matrix(Matrix::bdiag(dum_pen, pen_mat)) ) # add 0 matrix for non-penalized non-smooths and intercept
}

# cross validation
fun.gee1step.cv.exact <- function(w, d, grid, data,
                            namesd,
                            cname_,
                            fit.initial,
                            cv = FALSE, # cluster hold-1-out CV
                            K = 10,
                            B = 1000,
                            seed = 1){

  beta <- fit.initial$coefficients
  clust.vec <- as.vector(data[,..cname_])[[1]]
  clusters <- unique(clust.vec)
  n <- length(clusters)
  K <- min(c(K,n)) # make sure not too many folds
  N <- nrow(data)
  yindex <- unique(fit.initial$model$yindex.vec)
  L <- length(yindex)
  mse.ls <- NULL
  ###
  if(cv == TRUE | cv == "fastCV"){
    clust.ls <- lapply(clusters, function(ii) which(clust.vec == ii))
    YY <- as.matrix(data[, c("Y"), with = FALSE]) # outcome
    XX <- as.matrix(data[, ..namesd, with = FALSE]) # remove outcome

    # transform into list of held-out-sets
    YY <- lapply(clusters, function(ii) YY[ clust.ls[[ii]], ])
    XX <- lapply(clusters, function(ii) XX[ clust.ls[[ii]], ])
    ########################
    # hold-1-cluster out CV
    ########################
    if(is.list(grid)){
      grid.len <- length(grid)
      grid.ls <- grid
    }else{
      grid.ls <- list(grid)
      grid.len <- 1
    }

    for(ge in 1:grid.len){
      # iterate over grids
      grid <- grid.ls[[ge]]
      grid.size <- nrow(grid)
      mse.mat <- matrix(NA, nrow = nrow(grid), ncol = n)

      for(clust in 1:n){
        message(paste("cluster cv:", clust, "of ", n))

        di <- Reduce("+", d[-clust])  # sum over sample of clusters
        wi <- Reduce("+", w[-clust]) # sum over sample of clusters
        idx <- clust.ls[[clust]]
        n_i <- length(idx)
        fit.penal <- NULL

        for(g in 1:nrow(grid)){
          # construct penalty matrix
          penalty <- penalty.construct(model=fit.initial,
                                       lambda = grid[g,] * (N-n_i)/N, # adjust for sample size differences
                                       nonSmooths = 0) # nonSmooths set to 0 here

          # penalized estimator
          beta.boot.hat <- sanic::solve_chol(a=wi+penalty, b=di)
          mse.mat[g, clust] <- mean( (YY[[clust]] - XX[[clust]] %*% beta.boot.hat)^2 ) # fitted values on full range of covariates

        }
      }

      bias <- asy.var <- rep(NA, nrow(grid))
      mse <- as.numeric(rowMeans(mse.mat))

      if(ge < grid.len){
        # update new grid
        lambda.star <- grid[which.min(mse),] # optimal value
        grid.ls[[ge+1]] <- expand.grid(lapply(lambda.star, function(x) x*grid.ls[[ge+1]]))
      }
    }

  }else if(cv == "kfold"){
    YY <- as.matrix(data[, c("Y"), with = FALSE]) # outcome
    XX <- as.matrix(data[, ..namesd, with = FALSE]) # remove outcome

    # create folds
    indices <- 1:n
    set.seed(seed)
    sets <- split(indices[sample.int(n)], sort(rank(indices) %% K))
    clust.ls <- lapply(clusters, function(ii) which(clust.vec == ii))
    clust.ls <- lapply(1:K, function(kk) do.call(c,
                                                 clust.ls[ sets[[kk]] ] ) )

    # transform into list of held-out-sets
    YY <- lapply(1:K, function(ii) YY[ clust.ls[[ii]], ])
    XX <- lapply(1:K, function(ii) XX[ clust.ls[[ii]], ])

    # precompute fold specific quantities
    d.ls <- lapply(1:K, function(kk) Reduce("+", d[ -sets[[kk]] ]) ) #lapply(d, function(x) Reduce("+", x ) ) # sum over sample of clusters
    w.ls <- lapply(1:K, function(kk) Reduce("+", w[ -sets[[kk]] ]) ) #lapply(w, function(x) Reduce("+", x[-clust] ) ) # sum over sample of clusters

    ########################
    # k-fold CV
    ########################
    if(is.list(grid)){
      grid.len <- length(grid)
      grid.ls <- grid
    }else{
      grid.ls <- list(grid)
      grid.len <- 1
    }

    for(ge in 1:grid.len){
      # iterate over grids
      grid <- grid.ls[[ge]]
      grid.size <- nrow(grid)
      mse.mat <- matrix(NA, nrow = nrow(grid), ncol = K)
      for(kk in 1:K){
        message(paste("k-fold cv:", kk, "of ", K))
        idx <- clust.ls[[kk]] # rows of held out set
        n_i <- length(idx)
        fit.penal <- NULL
        wi <- w.ls[[kk]]
        di <- d.ls[[kk]]
        for(g in 1:nrow(grid)){
          # construct penalty matrix
          penalty <- penalty.construct(model=fit.initial,
                                       lambda = grid[g,] * (N-n_i)/N, # adjust for sample size differences
                                       nonSmooths = 0) # nonSmooths set to 0 here

          # penalized estimator
          beta.boot.hat <- sanic::solve_chol(a=wi+penalty, b=di) # coefficient estimate
          mse.mat[g, kk] <- mean( (YY[[kk]] - XX[[kk]] %*% beta.boot.hat)^2 ) # fitted values on full range of covariates
        }
      }

      bias <- asy.var <- rep(NA, nrow(grid))
      mse <- as.numeric(rowMeans(mse.mat))

      if(ge < grid.len){
        # update new grid
        lambda.star <- grid[which.min(mse),] # optimal value
        grid.ls[[ge+1]] <- expand.grid(lapply(lambda.star, function(x) x*grid.ls[[ge+1]]))
      }
    }
  }else if(cv == "fastkfold"){
    message("fast kfold exact")
    YY <- as.matrix(data[, c("Y"), with = FALSE]) # outcome
    XX <- as.matrix(data[, ..namesd, with = FALSE]) # remove outcome
    fit.penal <- NULL

    # create folds
    indices <- 1:n
    # rows of held out set
    set.seed(seed)
    sets <- split(indices[sample.int(n)], sort(rank(indices) %% K))
    clust.ls <- lapply(clusters, function(ii) which(clust.vec == ii))
    clust.ls <- lapply(1:K, function(kk) do.call(c,
                                                 clust.ls[ sets[[kk]] ] ) )

    # transform into list of held-out-sets
    YY <- lapply(1:K, function(ii) YY[ clust.ls[[ii]], ])
    XX <- lapply(1:K, function(ii) XX[ clust.ls[[ii]], ])

    ho.len <- sapply(sets, length)
    samp.vec <- sapply(clust.ls, length)

    # precompute fold specific quantities
    D <- lapply(1:K, function(kk) Reduce("+", d[ -sets[[kk]] ]) ) #lapply(d, function(x) Reduce("+", x ) ) # sum over sample of clusters
    D <- lapply(1:K, function(kk) D[[kk]] * N / (N-samp.vec[kk]) ) # scale D by sample size adjustment
    D <- do.call(cbind, D)
    wi <- Reduce("+", w) #* samp.scl

    ########################
    # fast k-fold CV
    ########################
    if(is.list(grid)){
      grid.len <- length(grid)
      grid.ls <- grid
    }else{
      grid.ls <- list(grid)
      grid.len <- 1
    }

    mse.ls <- vector(length = grid.len, "list")
    for(ge in 1:grid.len){
      # iterate over grids
      grid <- grid.ls[[ge]]
      grid.size <- nrow(grid)
      mse.mat <- matrix(NA, nrow = nrow(grid), ncol = K)
      for(g in 1:nrow(grid)){
        message(paste("k-fold cv:", g, "of ", nrow(grid)))
        # construct penalty matrix
        penalty <- penalty.construct(model=fit.initial,
                                     lambda = grid[g,], #  * (N-n_i)/N # adjust for sample size differences
                                     nonSmooths = 0) # nonSmooths set to 0 here

        # betas for all folds
        beta.boot.hat <- sanic::solve_chol(a = wi + penalty, b = D)

        for(kk in 1:K){
          mse.mat[g, kk] <- crossprod(YY[[kk]] - XX[[kk]] %*% beta.boot.hat[,kk]) / samp.vec[kk] # fitted values on full range of covariates
        }
      }

      bias <- asy.var <- rep(NA, nrow(grid))
      mse <- as.numeric(rowMeans(mse.mat))
      mse.ls[[ge]] <- cbind(grid, mse)

      if(ge < grid.len){
        # update new grid
        lambda.star <- grid[which.min(mse),] # optimal value
        grid.ls[[ge+1]] <- expand.grid(lapply(lambda.star, function(x) x*grid.ls[[ge+1]]))
      }
    }
  }

  lambda.star <- grid[which.min(mse),] # optimal value

  return(list(mse = mse,
              asy.var = asy.var,
              bias = bias,
              X = XX,
              grid = grid,
              lambda.star = lambda.star,
              fit.penal = fit.penal,
              mse.ls = mse.ls))

}

# cross validation
fun.gee1step.cv <- function(w, d, grid, data,
                            namesd,
                            cname_,
                            fit.initial,
                            cv = FALSE, # cluster hold-1-out CV
                            B = 1000,
                            K = 10,
                            seed = 1){

  data <- data.table::copy(data)
  link.fn <- fit.initial$family$link
  family.dist <- fit.initial$family$family

  # predictions
  if(link.fn == "identity"){
    pred.fn <- function(x){ x }
  }else if(link.fn == "logit"){
    pred.fn <- function(x){ 1 / (1 + exp(-x)) }
  }else if(link.fn == "log"){
    pred.fn <- function(x){ exp(x) }
  }

  # evaluations
  # https://stats.stackexchange.com/questions/71720/error-metrics-for-cross-validating-poisson-models
  # might consider the Dawid Sebastiani score -- seems fast an not super likelihood based
  if(family.dist == "gaussian"){
    criteria.fn <- function(preds, y){  crossprod(y-preds)/length(y) }# # MSE
  }else if(family.dist == "binomial"){
    criteria.fn <- function(preds, y){ -mean(dbinom(x=y, size=1, prob=preds, log=TRUE)) }
  }else if(family.dist == "poisson"){
    criteria.fn <- function(preds, y){ -mean(dpois(x=y, lambda=preds, log=TRUE)) }
  }

  beta <- fit.initial$coefficients
  # p <- length(beta)
  clust.vec <- as.vector(data[,..cname_])[[1]]
  clusters <- unique(clust.vec)

  n <- length(clusters)
  K <- min(c(n, K))
  N <- nrow(data)
  yindex <- unique(fit.initial$model$yindex.vec)
  L <- length(yindex)
  if(cv == TRUE){
    clust.ls <- lapply(clusters, function(ii) which(clust.vec == ii))
    YY <- as.matrix(data[, c("Y"), with = FALSE]) # outcome
    XX <- as.matrix(data[, ..namesd, with = FALSE]) # remove outcome

    # turn into list of held-out-sets
    YY <- lapply(clusters, function(ii) YY[ clust.ls[[ii]], ])
    XX <- lapply(clusters, function(ii) XX[ clust.ls[[ii]], ])
    fit.penal <- NULL

    ########################
    # hold-1-cluster out CV
    ########################
    if(is.list(grid)){
      grid.len <- length(grid)
      grid.ls <- grid
    }else{
      grid.ls <- list(grid)
      grid.len <- 1
    }

    for(ge in 1:grid.len){
      # iterate over grids
      grid <- grid.ls[[ge]]
      grid.size <- nrow(grid)
      mse.mat <- matrix(NA, nrow = nrow(grid), ncol = n)
      for(clust in 1:n){
        message(paste("cluster cv:", clust, "of ", n))
        di <- Reduce("+", d[-clust])  # sum over sample of clusters
        wi <- Reduce("+", w[-clust])  # sum over sample of clusters
        idx <- clust.ls[[clust]] # rows of held out set
        n_i <- length(idx)
        for(g in 1:nrow(grid)){
          # construct penalty matrix
          penalty <- penalty.construct(model=fit.initial,
                                       lambda = grid[g,] * (N-n_i)/N, # adjust for sample size differences
                                       nonSmooths = 0) # nonSmooths set to 0 here

          # penalized estimator
          inf_fn <- sanic::solve_chol(a=wi/(n-1) + penalty,
                                      b=di - penalty %*% beta * (n-1))  # coefficient estimate
          beta.boot.hat <- beta + inf_fn / (n-1) # coefficient estimate
          mse.mat[g, clust] <- criteria.fn(preds = pred.fn(XX[[clust]] %*% beta.boot.hat),
                                           y = YY[[clust]]) # fitted values on full range of covariates
          rm(inf_fn, beta.boot.hat)
          }
      }

      bias <- asy.var <- rep(NA, nrow(grid))
      mse <- as.numeric(rowMeans(mse.mat))
    }

    if(ge < grid.len){
      # update new grid
      lambda.star <- grid[which.min(mse),] # optimal value
      grid.ls[[ge+1]] <- expand.grid(lapply(lambda.star, function(x) x*grid.ls[[ge+1]]))
    }

  }else if(cv == "fastCV"){
    clust.ls <- lapply(clusters, function(ii) which(clust.vec == ii))
    samp.vec <- sapply(clust.ls, length)
    YY <- as.matrix(data[, c("Y"), with = FALSE]) # outcome
    XX <- as.matrix(data[, ..namesd, with = FALSE]) # remove outcome

    # turn into list of held-out-sets
    YY <- lapply(clusters, function(ii) YY[ clust.ls[[ii]], ])
    XX <- lapply(clusters, function(ii) XX[ clust.ls[[ii]], ])
    fit.penal <- NULL
    w <- Reduce("+", w)

    ########################
    # hold-1-cluster out CV - fast
    ########################
    if(is.list(grid)){
      grid.len <- length(grid)
      grid.ls <- grid
    }else{
      grid.ls <- list(grid)
      grid.len <- 1
    }

    for(ge in 1:grid.len){
      # iterate over grids
      grid <- grid.ls[[ge]]
      grid.size <- nrow(grid)
      mse.mat <- matrix(NA, nrow = nrow(grid), ncol = n)
      for(g in 1:grid.size){
        message(paste("cluster FASTcv:", g, "of ", grid.size))

        # construct penalty matrix
        penalty <- penalty.construct(model=fit.initial,
                                     lambda = grid[g,],
                                     nonSmooths = 0) # nonSmooths set to 0 here
        penalty.vec <- penalty %*% beta * n # penalty vector

        # penalized estimator
        wwi <- sanic::solve_chol( w/n + penalty )
        for(clust in 1:n){
          # in future, could try correction (2/23/25): n/(n-1)
          di <- Reduce("+", d[-clust]) * N / sum(samp.vec[-clust])  #added in n.bar to deal with different sample sizes, n_i (2/1/25)
          inf_fn <- wwi %*% (di - penalty.vec)  # coefficient estimate
          beta.boot.hat <- beta + inf_fn / n # (n-1) # coefficient estimate
          mse.mat[g, clust] <- criteria.fn(preds = pred.fn(XX[[clust]] %*% beta.boot.hat), # fitted values on full range of covariates
                                           y = YY[[clust]])
        }
      }

      bias <- asy.var <- rep(NA, nrow(grid))
      mse <- as.numeric(rowMeans(mse.mat))
    }

    if(ge < grid.len){
      # update new grid
      lambda.star <- grid[which.min(mse),] # optimal value
      grid.ls[[ge+1]] <- expand.grid(lapply(lambda.star, function(x) x*grid.ls[[ge+1]]))
    }

  }else if(cv == "kfold"){
    YY <- as.matrix(data[, c("Y"), with = FALSE]) # outcome
    XX <- as.matrix(data[, ..namesd, with = FALSE]) # remove outcome

    # create folds
    K <- ifelse(K > n, n, K) # make sure K is not larger than n
    indices <- 1:n
    sets <- split(indices[sample.int(n)], sort(rank(indices) %% K))
    clust.ls <- lapply(clusters, function(ii) which(clust.vec == ii))
    clust.ls <- lapply(1:K, function(kk) do.call(c,
                                clust.ls[ sets[[kk]] ] ) )
    ho.len <- sapply(sets, length)

    # precompute fold specific quantities
    d.ls <- lapply(1:K, function(kk) Reduce("+", d[ -sets[[kk]] ]) ) #lapply(d, function(x) Reduce("+", x ) ) # sum over sample of clusters
    w.ls <- lapply(1:K, function(kk) Reduce("+", w[ -sets[[kk]] ]) ) #lapply(w, function(x) Reduce("+", x[-clust] ) ) # sum over sample of clusters

    # turn into list of held-out-sets
    YY <- lapply(1:K, function(ii) YY[ clust.ls[[ii]], ])
    XX <- lapply(1:K, function(ii) XX[ clust.ls[[ii]], ])
    ########################
    # k-fold CV
    ########################
    if(is.list(grid)){
      grid.len <- length(grid)
      grid.ls <- grid
    }else{
      grid.ls <- list(grid)
      grid.len <- 1
    }

    for(ge in 1:grid.len){
      message(paste("fast k-fold cv list:", ge, "of ", grid.len))
      # iterate over grids
      grid <- grid.ls[[ge]]
      grid.size <- nrow(grid)
      mse.mat <- matrix(NA, nrow = nrow(grid), ncol = K)
      for(kk in 1:K){
        idx <- clust.ls[[kk]] # rows of held out set
        n_i <- length(idx)
        fit.penal <- NULL
        n.up <- n - ho.len[[kk]]
        wi <- w.ls[[kk]]
        di <- d.ls[[kk]]
        for(g in 1:nrow(grid)){
          # construct penalty matrix
          penalty <- penalty.construct(model=fit.initial,
                                       lambda = grid[g,] * (N-n_i)/N, # adjust for sample size differences
                                       nonSmooths = 0) # nonSmooths set to 0 here

          # penalized estimator
          wwi <- wi/n.up + penalty
          ddi <- di - penalty %*% beta * n.up
          inf_fn <- sanic::solve_chol(a=wwi, b=ddi)
          beta.boot.hat <- beta + inf_fn / n.up # coefficient estimate
          mse.mat[g, kk] <- criteria.fn(preds = pred.fn(XX[[kk]] %*% beta.boot.hat), # fitted values on full range of covariates
                                        y = YY[[kk]])
        }
      }

      bias <- asy.var <- rep(NA, nrow(grid))
      mse <- as.numeric(rowMeans(mse.mat))

      if(ge < grid.len){
        # update new grid
        lambda.star <- grid[which.min(mse),] # optimal value
        grid.ls[[ge+1]] <- expand.grid(lapply(lambda.star, function(x) x*grid.ls[[ge+1]]))
      }
    }

  }else if(cv == "fastkfold"){
    # could still use separate k-specific initial estimates for each fold:
    # https://stackoverflow.com/questions/78331965/updating-a-fitted-r-mgcvbam-model-reports-an-invalid-type-closure-error
    YY <- as.matrix(data[, c("Y"), with = FALSE]) # outcome
    XX <- as.matrix(data[, ..namesd, with = FALSE]) # remove outcome

    # create folds
    indices <- 1:n
    set.seed(seed)
    sets <- split(indices[sample.int(n)], sort(rank(indices) %% K))
    clust.ls <- lapply(clusters, function(ii) which(clust.vec == ii))
    clust.ls <- lapply(1:K, function(kk) do.call(c,
                                                 clust.ls[ sets[[kk]] ] ) )

    # turn into list of held-out-sets
    YY <- lapply(1:K, function(ii) YY[ clust.ls[[ii]], ])
    XX <- lapply(1:K, function(ii) XX[ clust.ls[[ii]], ])

    # ho.len <- sapply(sets, length)
    samp.vec <- sapply(clust.ls, length)

    # precompute fold specific quantities
    D <- lapply(1:K, function(kk) Reduce("+", d[ -sets[[kk]] ]) ) # sum over sample of clusters
    D <- lapply(1:K, function(kk) D[[kk]] * N / (N-samp.vec[kk]) ) # scale D by sample size adjustment
    D <- as.matrix(do.call(cbind, D))
    wi <- Reduce("+", w) / n
    fit.penal <- NULL

    ########################
    # fast k-fold CV
    ########################
    if(is.list(grid)){
      grid.len <- length(grid)
      grid.ls <- grid
    }else{
      grid.ls <- list(grid)
      grid.len <- 1
    }

    for(ge in 1:grid.len){
      message(paste("fast k-fold cv list:", ge, "of ", grid.len))
      # iterate over grids
      grid <- grid.ls[[ge]]
      grid.size <- nrow(grid)
      mse.mat <- matrix(NA, nrow = grid.size, ncol = K)

      for(g in 1:grid.size){
        # construct penalty matrix
        penalty <- penalty.construct(model=fit.initial,
                                     lambda = grid[g,], #  * (N-n_i)/N # adjust for sample size differences
                                     nonSmooths = 0) # nonSmooths set to 0 here

        # penalized estimator
        penalty_vec <- as.numeric(penalty %*% beta * n)
        # calculate influence function for all folds
        inf_fn.mat <- sanic::solve_chol(a = wi + penalty,
                                        b = D - penalty_vec)

        for(kk in 1:K){
          beta.boot.hat <- beta + inf_fn.mat[,kk] / n # coefficient estimate
          mse.mat[g, kk] <- criteria.fn(preds = pred.fn(XX[[kk]] %*% beta.boot.hat), # fitted values on full range of covariates
                                        y = YY[[kk]])
        }
        rm(inf_fn.mat, beta.boot.hat)
      }

      bias <- asy.var <- rep(NA, nrow(grid))
      mse <- as.numeric(rowMeans(mse.mat))

      if(ge < grid.len){
        # update new grid
        lambda.star <- grid[which.min(mse),] # optimal value
        # **** Consider adding ifelse(lambda.star == 0, 1, lambda.star) to prevent setting to 0
        grid.ls[[ge+1]] <- expand.grid(lapply(lambda.star, function(x) x*grid.ls[[ge+1]]))
      }
    }

  }

  lambda.star <- grid[which.min(mse),] # optimal value

  return(list(mse = mse,
              asy.var = asy.var,
              bias = bias,
              X = XX,
              grid = grid,
              lambda.star = lambda.star
              ))

}


# calculate rho for exchangeable correlation structure
exchangeable.rho <- function(dx){
  ### Estimate ICC

  drho <- dx[,
             list(
               .N,
               sum_r = sum( resid ),
               uss_r = sum( resid^2 )
             ), keyby = .(cname_, yindex.vec)]

  drho[, wt_ij := ( N * (N-1) / 2), by = yindex.vec]
  drho[, rho_ij := sum_r^2 - uss_r, by = yindex.vec ]
  rho <- drho[, (sum(rho_ij)/2) / sum(wt_ij), by = yindex.vec] # GL: not sure why there is a division by two here
  rho_vec <- as.numeric(rho$V1[order(rho$yindex.vec)])
  rho_vec <- ifelse(rho_vec > 1, 0.99, rho_vec)
  rho_vec <- ifelse(rho_vec < -1, -0.99, rho_vec)

  # join
  colnames(rho)[2] <- "rho"
  return(dx[rho, on = .(yindex.vec)]) # add rho
}

