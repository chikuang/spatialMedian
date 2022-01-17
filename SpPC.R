# source("gsj.R")
# %function [pc,wpc,lambda,scores] = SpPC(t,x,q,wmu,dtyp)
# %
# %Functional Spherical Principal Components
# %
#' @param
# %INPUT: See comments on function SpMed.
# % Input variables:
#   %   t: Time grid (1 x m vector). Must be SORTED (in increasing order)
# %   x: Discretized curves (as in SpMed).
# %   q: Number of PCs to estimate (scalar).
# %   wmu: Spatial median weights. This is output 'w' of SpMed (n x 1 vector).
# %   dtyp: Data type (as in SpMed; optional).
# %
#'
#' @output
#' pc: Spherical PCs (q x m matrix).
#' wpc: Weights (n x q matrix). (See paper for explanation).
#' lambda: Eigenvalues (q x 1 vector). (Note that these eigenvalues are NOT consistent estimators of the true model eigenvalues; use e.g. median of squared PC scores to estimate eigenvalues consistently. See paper for further details).
# %   scores: Individual PC scores (n x q matrix). Useful for computing
# %   consistent estimators of the eigenvalues and for curve reconstruction
# %   (x ~ ones(n,1)*med + scores*pc)
# %
# % External functions called: gsj.m
# %

SpPC <- function(tt, x, q, wmu, dtyp = "s"){
  tt <- as.matrix(tt, ncol = length(tt), nrow = 1)
  eps <- 2^(-52)
  if (nargs() < 4)
    stop("Not enough  input variables")

  n <- nrow(x)
  m <- ncol(x)
  if (nrow(tt) > 1){
    tt <- t(tt)
  }
  if (ncol(tt) != m)
    stop("Dimensions of T and X not compatible")
    
  if(nrow(wmu) != n){
    wmu <- t(wmu)
  }
  if(nrow(wmu) != n)
    stop("Dimensions of WMU and X not compatible")
  
  if (!(tolower(dtyp) != 'n' | tolower(dtyp) != 's'))
    stop('Wrong input value for DTYP')
  
  # Inner-product matrices
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
  
  B <- (diag(1, nrow = n, ncol = n) - matrix(1, nrow = n, ncol = 1) %*% t(wmu)) %*% A %*%
    (diag(1, nrow = n, ncol = n) - wmu %*% matrix(1, nrow = 1, ncol = n))
  norms <- sqrt(diag(B))
  DN <- matrix(0, nrow= n, ncol = 1)
  I <- which(norms > eps)
  DN[I] <- 1/norms[I]
  B <- diag(as.numeric(DN)) %*% B %*% diag(as.numeric(DN))
  
  B <- (B+t(B))/2 # just for the numerical reasons
  
  
  # PC computation
  # temp <- eigen(B)
  # V <- temp$vectors[,1:q]
  # 
  # D <- temp$values[1:q]
  # 
  
  temp <- RSpectra::eigs_sym(B, q, which = "LA")
  V <- temp$vectors
  D <- temp$value
  # Check if the -V is greater
  
  
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
  }
  
  pc <- t(wpc) %*% diag(as.numeric(DN)) %*% (diag(1, n) - matrix(1, nrow = n, ncol = 1) %*% t(wmu)) %*% x
  scores <- 0.5 * ( pc[, 1:(m-1)] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% 
                      t((diag(1, n) - matrix(1, n, 1) %*% t(wmu)) %*% x[, 1:(m-1)]) +
                      pc[, 2:m] %*% diag(tt[2:m] - tt[1:(m-1)]) %*% 
                      t((diag(1, n) - matrix(1, n, 1) %*% t(wmu)) %*% x[, 2:m]))
  
  scores <- t(scores)

  
  return(list(pc = pc, wpc = wpc, lambda = lambda, scores = scores))
}