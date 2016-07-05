#' thetaltivariate Laplace sithetalation function.
#'
#' Function that sithetalates from the thetaltivariate Laplace distributions.
#'
#' @param m Numeric value, length of t vector.
#' @param theta Numeric vector with length m.
#' @param theta Numeric vector with length m.
#' @param kappa Numeric value with length m.
#' @param sigma Numeric value (positive).
#' @param Sigma Matrix (m times m).
#' @param tau Numeric value (positive).
rMultLaplace = function(m, theta, gamma, Sigma, sigma = 0.2, tau = 1) {

  if (length(theta) != m) stop(paste('theta must have length', m))
  if (length(gamma) != m) stop(paste('gamma must have length', m))
  if (sigma < 0) stop('sigma must be positive')
  if (!is.matrix(Sigma))
    stop('Sigma must be matrix')
  else
    if (nrow(Sigma) != m | ncol(Sigma) != m) stop(paste('Sigma must have dimension', m, "times", m))

#  m  = length(t)
  x0 = rnorm(m, sd = sigma)
  w  = rgamma(1, shape = tau)

  x   = t(chol(Sigma)) %*% x0
  out = as.numeric(theta + sqrt(w) * x + w * gamma)

  attr(out, "x") = as.numeric(x)
  attr(out, "w") = as.numeric(w)

  return(out)
}




