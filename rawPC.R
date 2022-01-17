# source("gsj.R")

#' Functional Principal Components (no roughness penalization)
#' 
#' @parameter
#' 
#' t: Time grid (1 x m vector). Must be SORTED (in increasing order).
#' x: Discretized curves (n x m matrix). Each row represents a different
#'   curves. Missing values NOT allowed. If your data set has missing
#'   values, use either interpolation or smoothing to fill in the gaps.
#' q: Number of PCs to estimate (scalar).
#' dtyp: Data type (optional; character variable). Enter 's' if the data 
#'         is "reasonably" smooth, or 'n' if the data is very noisy (in that case a bias correction is necessary and the error variance is estimated nonparametrically using the GSJ estimator; this will slow things down, so avoid if possible. It may also produce an inner product matrix that is non-negative definite, in which case the uncorrected matrix will be used).  Default value: 's'.
#'         
#' @return 
#' pc: Raw PCs (q x m matrix)
#' lambda: Eigenvalues (q x 1 vector)
#' scores: Individual PC scores (n x q matrix). Useful for individual curve estimation (x ~ ones(n,1)*mean(x) + scores*pc)
#' 
#' @note external function "gsj.R"

RawPC <- function(tt, x, q, dtyp = "s"){
  eps <- 2^(-52)
  # Input check
  
  if (nargs() < 3)
    stop("Not enough  input variables")
  n <- nrow(x)
  m <- ncol(x)

  if(length(tt) != m)
    stop('Dimensions of T and X not compatible')
  
  if (!(tolower(dtyp) != 'n' | tolower(dtyp) != 's'))
    stop('Wrong input value for DTYP')
  
  # Inner-product matrices
  A <- 0.5 * (x[, 1:(m-1)] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% t(x[, 1:(m-1)]) +
                x[, 2:m] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% t(x[,2:m]))
  
  
  if( tolower(dtyp) == "n"){
    se2 <- gsj(kronecker(matrix(1, 1, n), tt), t(x))
    A1 <- A - diag(se2) * (tt[m] - tt[1]);
    if (min(eig(A1)>-1e-6))
      A <- A1
    else{
      print('Corrected inner-product matrix is not nonnegative definite \n Using uncorrected matrix instead ')
    }
  }
  B <- (diag(1, n) - matrix(1, n, n)/n) %*% A %*% (diag(1, n) - matrix(1, n, n)/n)
  B <- (B+t(B))/2 # just for the numerical reasons
  
  ## PC computation
  temp <- RSpectra::eigs_sym(B, q, which = "LA")
  V <- temp$vectors
  D <- temp$value
  
  if (is.ordered(D[seq(q, 1, by = -1)])){
    wpc <- V %*% diag(1/sqrt(D))
    lambda <- D/n
  } else if (is.ordered(D)) {
    wpc <- V[, seq(q, 1, by = -1)] %*% diag(1/sqrt(D[seq(q, 1, by = -1)]))
    lambda <-  D[seq(q, 1, by = -1)]/n
  } else{
    I <- order(D, decreasing = TRUE)
    D <- sort(D, decreasing = TRUE)
    
    wpc <- V[, I] %*% diag(1/sqrt(D))
    lambda <- D/n
    
    pc <- t(wpc) %*% (diag(1, n) - matrix(1, n, n)/n) %*% x
    scores <- 0.5 * (pc[,1:(m-1)] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% 
                       t((diag(1, n) - matrix(1, n, n)/n) %*% x[, 1:(m-1)] ) +
                       pc[, 2:m] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% 
                       t((diag(1, n) - matrix(1, n, n)/n) %*% x[, 2:m]))
    scores <- t(scores)
    return(list(pc = pc, lambda = lambda, scores = scores))
  }
}

# x <-  matrix(   c(0.1736,    0.9892,    0.2928,    0.3668,    0.1895,    0.8960,
# 0.5752,    0.4899,    0.8014,    0.7395,    0.1237,    0.5154,
# 0.6062,    0.6949,    0.3465,    0.5247,    0.8210,    0.5445,
# 0.2144,    0.4114,    0.0833,    0.8045,    0.6379,    0.6064,
# 0.5199,    0.0348,    0.5111,    0.8169,    0.0161,    0.7604), nrow = 5, ncol = 6, byrow = T)
# tt <- seq(0, 1, length.out = 6)
# q <- 4
# RawPC(tt, x, q)
