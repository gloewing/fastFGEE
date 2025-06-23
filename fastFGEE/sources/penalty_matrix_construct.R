pen_mat_construct <- function(model, nonSmooths = 1){

  # penalty
  S <- model$smooth
  L <- length(S)
  dum_pen <- matrix(0, ncol = nonSmooths+1, nrow = nonSmooths+1) # dummy terms
  # extract penalty matrices and scale by lambda -- verified (at least unweighted case)
  pen_mat <- lapply(1:L, function(i) S[[i]]$S[[1]] * model$sp[i] ) # * (model$sp[i]) / (S[[i]]$S.scale) #  * (model$sp[i]) / S[[i]]$S.scale # * sqrt(model$sp[i])
  pen_mat <- Matrix::bdiag(pen_mat) # concatenate matrices into block diagonal
  
  return( as.matrix(Matrix::bdiag(dum_pen, pen_mat)) ) # add 0 matrix for non-penalized non-smooths and intercept
}
 
# pen_mat_update <- function(model, nonSmooths = 1){
#   
#   # penalty
#   S <- model$smooth
#   L <- length(S)
#   dum_pen <- matrix(0, ncol = nonSmooths+1, nrow = nonSmooths+1) # dummy terms
#   # extract penalty matrices and scale by lambda -- verified (at least unweighted case)
#   pen_mat <- lapply(1:L, function(i) S[[i]]$S[[1]] * model$sp[i] ) # * (model$sp[i]) / (S[[i]]$S.scale) #  * (model$sp[i]) / S[[i]]$S.scale # * sqrt(model$sp[i])
#   pen_mat <- Matrix::bdiag(pen_mat) # concatenate matrices into block diagonal
#   
#   return( as.matrix(Matrix::bdiag(dum_pen, pen_mat)) ) # add 0 matrix for non-penalized non-smooths and intercept
# } 
  # smooth_dim <- sum(sapply(S, function(x) ncol(x$S[[1]]))) # number of dimensions of smooths
  # pen_dim <-  nonSmooths + smooth_dim
  

  # S <- model$smooth[[1]]$S[[1]] #sm[[1]]$S[[1]] ## penalty matrix
  # alpha <- sapply(model$smooth, "[[", "S.scale")
  # lambda <- alpha/model$sp  ## lambda parameter
  # smooth = model$smooth[[1]] # assumes only 1 penalty matrix
  # par = smooth$first.para:smooth$last.para
  # p <- length(beta_hat)
  # pen_mat <- matrix(0, ncol = p, nrow = p)
  # pen_mat[par,par] <- lambda * S # penalty matrix - absorb constant
  # 
  # 
  
