#' source("gsj.R")

#' %Functional Spatial Median
#' @note 
#' It is assumed that all curves are observed on the same time grid. The time grid can be irregular (doesn't need to be equispaced).
#' 
#' @note require external function "gsj.R"
#' 
#' @parameter
#' t: Time grid (1 x m vector). Must be SORTED (in increasing order).
#' x: Discretized curves (n x m matrix). Each row represents a different curve. Missing values NOT allowed. If your data set has missing values, use either interpolation or smoothing to fill in the gaps.
#' dtyp: Data type (optional; character variable). Enter 's' if the data is "reasonably" smooth, or 'n' if the data is very noisy (in that case a bias correction is necessary and the error variance is estimated nonparametrically using the GSJ estimator; this will slow things down, so avoid if possible. It may also produce an inner-product matrix that is not nonnegative definite, in which case the uncorrected matrix will be used). Default value: 's'.
#' @output
# %   med: Spatial median (1 x m vector).
# %   w: Weights (n x 1 vector). The spatial median is w'*x.

SpMed <- function(tt, x, dtyp = 's'){
  tt <- as.matrix(tt, ncol = length(tt), nrow = 1)
  eps <- 2^(-52)
  # Input check
  if (nargs() < 2)
    stop("Not enough input variables")
  
  n <- nrow(x)
  m <- ncol(x)

  if (nrow(tt) > 1)
    tt <- t(tt)
  
  if (ncol(tt) != m)
    stop('Dimensions of T and X not compatible')
  
  if (!(dtyp != 'n' | dtyp != 's'))
    stop('Wrong input value for DTYP')
  
  # Inner-product matrix
  A <- 0.5 * (x[, 1:(m-1)] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% t(x[, 1:(m-1)]) +
    x[, 2:m] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% t(x[, 2:m])) 
  
  if( tolower(dtyp) == "n"){
    se2 <- gsj(kronecker(matrix(1, 1, n), tt), t(x))
    A1 <- A - diag(se2) * (tt[m] - tt[1]);
    if (min(eig(A1)>-1e-6))
      A <- A1
    else{
      print('Corrected inner-product matrix is not nonnegative definite \n Using uncorrected matrix instead ')
    }
  }
  
  ## Iterative minimization from sample mean
  w <- matrix(1, nrow = n, ncol = 1)/n # a very naive way
  norms <- sqrt(diag(A) + as.vector(t(w) %*% A %*% w) - 2 * A %*% w) # pay attention to the R syntax that it does not like to do arithmetic on a number to vector
  f <- sum(norms)
  err <- 1
  iter <- 0
  
  while( err > 1E-5 && iter < 50){
    iter <- iter + 1
    f0 <- f
    if (min(norms < eps)){
      i0 = find(norms < eps)
      w <- matrix(0, nrow = n, ncol = 1)
      w[i0] <- 1/length(i0)
    } else{
      w <- 1/norms
      w <- w/sum(w)
    }
    
    norms <- sqrt(diag(A) + as.vector(t(w) %*% A %*% w) - 2 * A %*% w) # same here
    f <- sum(norms)
    err <- abs(f/f0 - 1)
  }
  
  med <- t(w) %*% x
  return(list(med = med, w = w, norms = norms))
}

# # test case
# xx <- matrix(c(    0.6896,    0.5850,    0.4926,    0.9760,    0.3091,    0.8975,    0.0293,    0.8467,    0.0540,    0.7962,0.7889,
#                    0.1318,    0.0734,    0.6549,    0.0364,    0.1209,    0.4996,    0.5279,    0.2461,    0.0206,    0.6179,0.0924,
#                    0.1235,    0.8223,    0.8901,    0.3262,    0.9158,    0.6153,    0.0321,    0.5815,    0.6815,    0.0702,0.2379,
#                    0.1909,    0.7229,    0.5385,    0.9730,    0.1355,    0.5831,    0.8271,    0.9377,    0.5986,    0.0693,0.2436,
#                    0.1457,    0.9259,    0.2822,    0.3650,    0.3321,    0.6983,    0.3400,    0.0478,    0.1140,    0.1360,0.1048), nrow = 5, byrow = T)
# 
# tt <- matrix(seq(0, 1, length.out = 11), ncol = 1)
# 
# ans <- SpMed(tt, xx)
# ww <- round(ans$w, 4)
# medd <- round(ans$med, 4)
