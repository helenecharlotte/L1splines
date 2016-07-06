#' @title Laplace cumulative distribution function.
#'
#' The cumulative distribution function for the univariate Laplace distributions.
#' The parametrization with kappa is used.
#'
#' @param x Numeric value, where the function is evaluated.
#' @param mu Numeric value. If specified, this parametrization will be used.
#' @param theta Numeric value.
#' @param kappa Numeric value (positive).
#' @param sigma Numeric value (positive).
pLaplace = function(x, mu=Inf, theta=0, kappa=1, sigma=1) {
  if (sigma < 0) stop("sigma must be positive")
  if (kappa < 0) stop("kappa must be positive")
  if (!(mu == Inf)) kappa = sqrt(2) * sigma / (mu + sqrt(2 * sigma^2 + mu^2))
  if (x < theta)
    out = kappa^2 / (1 + kappa^2) * exp(- sqrt(2) / (sigma * kappa) * abs(x - theta))
  else
    out = 1 - 1 / (1 + kappa^2) * exp(- sqrt(2) / (sigma * kappa) * abs(x - theta))
  return(out)
}

#' @title Laplace quantile function.
#'
#' The quantile function for the univariate Laplace distributions.
#' Obtained from numerically inverting pLaplace.
#' The parametrization with kappa is used
#'
#' @param q Numeric value, where the function is evaluated, between 0 and 1.
#' @param mu Numeric value. If specified, this parametrization will be used.
#' @param theta Numeric value.
#' @param kappa Numeric value (positive).
#' @param sigma Numeric value (positive).
qLaplace = function(q, mu=Inf, theta=0, kappa=1, sigma=1) {
  if (q < 0 || q > 1) stop("q must be in (0,1)")
  if (sigma < 0) stop("sigma must be positive")
  if (kappa < 0) stop("kappa must be positive")
  if (!(mu == Inf)) kappa = sqrt(2) * sigma / (mu + sqrt(2 * sigma^2 + mu^2))
  return(nleqslv(0, function(x)
    pLaplace(x, theta = theta, kappa = kappa, sigma = sigma) - q)$x)
}

#' @title Laplace density function.
#'
#' The density function for the univariate Laplace distributions.
#' Can use both the parametrization with respect to kappa and mu.
#'
#' @param x Numeric value, where the function is evaluated.
#' @param mu Numeric value. If specified, this parametrization will be used.
#' @param theta Numeric value.
#' @param kappa Numeric value (positive).
#' @param sigma Numeric value (positive).
dLaplace = function(x, mu=Inf, theta=0, kappa=1, sigma=1) {
  if (sigma < 0) stop("sigma must be positive")
  if (kappa < 0) stop("kappa must be positive")
  if (!(mu == Inf)) kappa = sqrt(2) * sigma / (mu + sqrt(2 * sigma^2 + mu^2))
  return(sqrt(2) / sigma * kappa / (1 + kappa^2) *
           exp(-sqrt(2) * kappa / sigma * abs(x - theta) * (x >= theta)) *
           exp(-sqrt(2) / (sigma * kappa) * abs(x - theta) * (x < theta)))
}

#' @title Laplace simulation function.
#'
#' Function that simulates from the univariate Laplace distributions.
#' Two parametrizations can be used.
#'
#' @param m Numeric value, length of t vector.
#' @param mu Numeric value. If specified, this parametrization will be used.
#' @param theta Numeric value.
#' @param kappa Numeric value (positive).
#' @param sigma Numeric value (positive).
rLaplace = function(m, mu=Inf, theta=0, kappa=1, sigma=1) {
  if (sigma < 0) stop("sigma must be positive")
  if (kappa < 0) stop("kappa must be positive")
  if (!(mu == Inf)) kappa = sqrt(2) * sigma / (mu + sqrt(2 * sigma^2 + mu^2))
  U1 = runif(m)
  U2 = runif(m)
  return(theta + sigma / sqrt(2) * log(U1^kappa / U2^(1 / kappa)))
}

