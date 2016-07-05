#' Split Bregman L1 spline.
#'
#' Function to compute the L1 spline, by use of the split Bregman method.
#'
#' @param x Numeric value(s) to evaluate the spline in.
#' @param t Numeric vector of measurement times.
#' @param y Numeric vector of values to be fitted.
#' @param gfun Choice of Greens function (choose from G1-G10).
#' @param lbratio Numeric value. Ratio of lambda over beta.
#' @param beta Numeric value.
#' @param niter Number of iterations.
#' @param stopiter Stopping criterion.
L1splineSB = function(x, t, y, gfun, lbratio=0.001, beta=1, niter=200, stopiter=1e-10) {

  m = length(t)

  d = z = b = matrix(0, m, niter)

  shrink = function(x, gamma)
    x / abs(x) * max(abs(x) - gamma, 0)

  dist = function(x, y)
    return(t(x-y) %*% (x-y))

  for (k in 1:(niter-1)) {

    ystar    = d[, k] + y - b[, k]
    z[, k+1] = thetafun(t, ystar, gfun, lambda=lbratio)
    d[, k+1] = sapply(z[, k+1] - y + b[, k], function(x) shrink(x, 1/(2*beta)))
    b[, k+1] = b[, k] + (z[, k+1] - y - d[, k+1])

    if (dist(z[, k+1],  z[, k]) < stopiter) break
  }

  ystar = d[, k] + y - b[, k]

  out = thetafun(t, ystar, gfun, x = x, lambda=lbratio)
  attr(out, "niter") = k

  return(out)
}

L1splineSB = function(x, t, y, gfun, lbratio=0.001, beta=1, niter=200, stopiter=1e-10) {

  m = length(t)

  d = z = b = matrix(0, m, niter)

  shrink = function(x, gamma)
    x / abs(x) * max(abs(x) - gamma, 0)

  dist = function(x, y)
    return(t(x-y) %*% (x-y))

  for (k in 1:(niter-1)) {

    ystar    = d[, k] + y - b[, k]
    z[, k+1] = thetafun(t, ystar, gfun, lambda=lbratio)
    d[, k+1] = sapply(z[, k+1] - y + b[, k], function(x) shrink(x, 1/(2*beta)))
    b[, k+1] = b[, k] + (z[, k+1] - y - d[, k+1])

    if (dist(z[, k+1],  z[, k]) < stopiter) break
  }

  ystar = d[, k] + y - b[, k]

  out = thetafun(t, ystar, gfun, x = x, lambda=lbratio)
  attr(out, "niter") = k

  return(out)
}

