#' Adapted L1 spline.
#'
#' Function to compute the adapted L1 spline, by use of the continuation method.
#'
#' @param x Numeric value(s) to evaluate the spline in.
#' @param t Numeric vector of measurement times.
#' @param y Numeric vector of values to be fitted.
#' @param gfun Choice of Greens function (choose from G1-G10).
#' @param sigma Variance, per default set to srqt(2). Effect translated into lambda.
#' @param b Truncation window.
#' @param lambda Smoothness parameter.
#' @param niter Number of iterations.
#' @param stopiter Stopping criterion.
L1splineAdapted = function(x, t, y, gfun, sigma=sqrt(2), b=0.5, lambda=0.0001, niter=200, stopiter=1e-10) {

  m    = length(t)
  beta = 1.1

  theta = psi = matrix(y, m, niter)

  L1opt = function(x) {
    out = y
    out[x > y + sqrt(2)/(2*sigma*beta)] = (-sqrt(2)/(2*sigma*beta) + x)[x > y + sqrt(2)/(2*sigma*beta)]
    out[x < y - sqrt(2)/(2*sigma*beta)] = (sqrt(2)/(2*sigma*beta) + x)[x < y - sqrt(2)/(2*sigma*beta)]
    return(out)
  }

  L1optTrunc = function(x) {
    out = y
    out[x > y + sqrt(2)/(2*sigma*beta)] = (-sqrt(2)/(2*sigma*beta) + x)[x > y + sqrt(2)/(2*sigma*beta)]
    out[x < y - sqrt(2)/(2*sigma*beta)] = (sqrt(2)/(2*sigma*beta) + x)[x < y - sqrt(2)/(2*sigma*beta)]
    out[abs(x - y) < b]                 = x[abs(x - y) < b]
    return(out)
  }

  dist = function(x, y)
    return(t(x-y) %*% (x-y))

  if (niter > 1) {
    for (k in 1:(niter-1)) {

      beta     = beta * 1.1

      theta[, k+1] = thetafun(t, psi[, k], gfun, lambda=lambda/beta)
      (psi[, k+1]  = L1optTrunc(theta[, k+1]))

      if (k > 1 & dist(theta[, k+1],  theta[, k]) < stopiter) break
    }
    yout = theta[, k+1]
  } else {
    yout = y
    k    = 1
  }

  out = thetafun(t, yout, gfun, x = x, lambda=lambda)

  attr(out, "niter") = k
  attr(out, "sigma") = sigma

  return(out)
}

