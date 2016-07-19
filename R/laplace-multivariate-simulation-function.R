#' @title Laplace simulation function.
#'
#' Function that simulates from the univariate Laplace distributions.
#' Two parametrizations can be used.
#'
#' @param n Numeric value, number of simulations.
#' @param theta Numeric vector of position parameters of length d.
#' @param gamma Numeric vector of skewness parameters of length d.
#' @param Sigma Numeric positive definite matrix of size d times d.
#' @param sigma Numeric value (positive).
#' @param tau Numeric value (positive).
#' @param seed Integer value to be specified if want certain seed.
#' @export
rMultiLaplace = function(n=1, theta=c(0, 0), gamma=c(0,0), Sigma=diag(2),
                         sigma=0.2, tau=1, seed=sample(2500, 1)) {
  if (sigma < 0) stop("sigma must be positive")
  if (tau < 0) stop("tau must be positive")
  if (!is.positive.definite(Sigma)) stop("Sigma must be positive definite")
  if (length(theta) != length(gamma) | length(theta) != dim(Sigma)[1] | length(gamma) != dim(Sigma)[1])
    stop("check dimensions of parameters")
  ##--- simulation of multivariate Gaussian
  rMultiGaussian = function(Sigma, seed=sample(2500, 1)) {
    set.seed(seed)
    m0 = dim(Sigma)[1]
    x0 = rnorm(m0, sd = sigma)
    return(t(chol(Sigma)) %*% x0)
  }

  ##--- simulation of multivariate Laplace
  out = sapply(1:n, function(i) {
    X  = rMultiGaussian(Sigma, seed=seed+i)
    W  = rgamma(1, shape = tau)
    Y  = theta + sqrt(W) * X + W * gamma
    return(Y)
  })

  return(out)
}
