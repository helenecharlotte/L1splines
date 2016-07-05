#' L1 spline (unknown dynamic variance).
#'
#' Function to compute an L1 spline, by use of the continuation method.
#' Fitting variance, allowing for dynamical change.
#'
#' @param x Numeric value(s) to evaluate the spline in.
#' @param t Numeric vector of measurement times.
#' @param y Numeric vector of values to be fitted.
#' @param gfun Choice of Greens function (choose from G1-G10).
#' @param t0 Anchor points for changing kappa.
#'        Per default set a division of [0,1] based on the number of observations.
#' @param sigma Variance, per default set to srqt(2). Effect translated into lambda.
#' @param div Integer determining where to put boundary knots.
#' @param lambda Smoothness parameter.
#' @param niter Number of iterations.
#' @param stopiter Stopping criterion.
L1splineSD = function(x, t, y, gfun, t0 = Inf, div = 10,
                      lambda=0.0001, niter=200, stopiter=1e-10, ...) {

  Qs = function(theta, sigma, y, m)
    return(- m / 2 * log(2) + sum(log(sigma)) + sum(sqrt(2) / sigma * abs(y - theta)))

  sigmat = function(t, tvec, sigmavec) {
    if (any(sigmavec < 0)) stop('Negative variance')
    return(approx(tvec, sigmavec, xout = t, rule = 2)$y)
  }

  m0  = floor(m / 50 + 9)
  tm0 = seq(0, 1, length = m0)

  if (any(t0 == Inf)) t0 = tm0[seq(2, m0-1, by = 2)]

  newfun = function(x) {
    sigma = sigmat(t, t0, abs(x))
    theta = L1spline(t, t, y, gfun, sigma = sigma, lambda = lambda, niter = niter, stopiter = stopiter)
    return(Qs(theta, sigma, y, m))
  }

  sigmavec = abs(optim(rep(1, length(t0)), newfun)$par)
  sigmahat = sigmat(t, t0,
                    sigmavec)

  out = L1spline(x, t, y, gfun, sigma = sigmahat,
                 lambda = lambda, niter = niter, stopiter = stopiter)
  attr(out, "sigmavec") = sigmavec
  attr(out, "tw")       = t0
  return(out)
}
