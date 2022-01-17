#' Gasser-Sroka-Jennen-Steinmetz nonparametric variance estimator
#' @parameters:
#'  x: Design points (n,m)
#'  y: Function observations (n,m)
#' @note  
#' Here n is the sample size and m is the number of functions or datasets
#'   

gsj <- function(x, y){
  n <- nrow(x)
  m <- ncol(x)
  if (n == 1){
    x <-  t(x)
    y <- t(y)
    n <- nrow(x)
    m <- ncol(x)
  }
  
  if (n < 3)
    return(NULL)
  ind <- apply(x, 2, order)
  x <- apply(x, 2, sort)
  
  # Also sort the y accordingly using the indicator above
  for (j in 1:m) {
    y[ ,j] <- y[ind[ ,j], j]
  }
  
  a <- (x[3:n, ] - x[2:(n-1), ]) / (x[3:n, ] - x[1:(n-2), ])
  b <- (x[2:(n-1), ] - x[1:(n-2), ])/(x[3:n, ]-x[1:(n-2), ])
  c <- a^2 + b^2 + 1
  e <- a * y[1:(n-2), ] + b * y[3:n, ] - y[2:(n-1),]
  
  estvar <- colMeans(e^2/c)
  return(estvar)
}
