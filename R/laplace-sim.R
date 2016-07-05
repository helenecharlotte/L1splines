#' Laplace simulation function.
#'
#' Function that simulates from the univariate Laplace distributions.
#' Two parametrizations can be used.
#'
#' @param m Numeric value, length of t vector.
#' @param mu Numeric value.
#' @param theta Numeric value.
#' @param kappa Numeric value.
#'              If not specified: Use one parametrization.
#'              If specified: Use other parametrization.
rLaplace = function(m, mu = 0, theta = 0, kappa = Inf, sigma = 1) {
## Recall: mu = 0 corresponds to symmetric, L(theta, sigma)
  if (any(kappa == Inf))
    kappa = sqrt(2) * sigma / (mu + sqrt(2 * sigma^2 + mu^2))
  U1 = runif(m)
  U2 = runif(m)
  return(theta + sigma / sqrt(2) * log(U1^kappa / U2^(1 / kappa)))
}

