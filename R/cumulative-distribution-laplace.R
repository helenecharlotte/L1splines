#' @title Laplace cumulative distribution function.
#'
#' The cumulative distribution function for the univariate Laplace distributions.
#' The parametrization with kappa is used
#'
#' @param x Numeric value, where the function is evaluated.
#' @param theta Numeric value.
#' @param kappa Numeric value (positive).
#' @param sigma Numeric value (positive).
pLaplace = function(x, theta = 0, kappa = 1, sigma = 1) {
  ## Recall: kappa = 1 corresponds to symmetric, L(theta, sigma)
  if (sigma < 0) stop("sigma must be positive")
  if (kappa < 0) stop("kappa must be positive")
  if (x < theta)
    out = kappa^2 / (1 + kappa^2) * exp(- sqrt(2) / (sigma * kappa) * abs(x - theta))
  else
    out = 1 - 1 / (1 + kappa^2) * exp(- sqrt(2) / (sigma * kappa) * abs(x - theta))
  return(out)
}

#' Laplace quantile function.
#'
#' The quantile function for the univariate Laplace distributions.
#' Obtained from numerically inverting pLaplace.
#' The parametrization with kappa is used
#'
#' @param q Numeric value, where the function is evaluated, between 0 and 1.
#' @param theta Numeric value.
#' @param kappa Numeric value (positive).
#' @param sigma Numeric value (positive).
qLaplace = function(q, theta = 0, kappa = 1, sigma = 1) {
  ## Recall: kappa = 1 corresponds to symmetric, L(theta, sigma)
  if (q < 0 || q > 1) stop("q must be in (0,1)")
  if (sigma < 0) stop("sigma must be positive")
  if (kappa < 0) stop("kappa must be positive")
  return(nleqslv(0, function(x)
    pLaplace(x, theta = theta, kappa = kappa, sigma = sigma) - q)$x)
}

#' Laplace density function.
#'
#' The density function for the univariate Laplace distributions.
#' Can use both the parametrization with respect to kappa and mu.
#'
#' @param x Numeric value, where the function is evaluated.
#' @param theta Numeric value.
#' @param mu Numeric value.
#' @param kappa Numeric value (positive).
#'        if not specified, parametrization with mu is used.
#'        if specified (numeric value), parametrization with kappa used.
#' @param sigma Numeric value (positive).
dLaplace = function(x, mu=1, sigma=0.5, kappa = Inf, theta=0.5) {
  if (kappa == Inf) kappa = sqrt(2) * sigma / (mu + sqrt(2 * sigma^2 + mu^2))
  return(sqrt(2) / sigma * kappa / (1 + kappa^2) *
           exp(-sqrt(2) * kappa / sigma * abs(x - theta) * (x >= theta)) *
           exp(-sqrt(2) / (sigma * kappa) * abs(x - theta) * (x < theta)))
}
