#' L1 spline.
#'
#' Function to compute the L1 spline, by use of the continuation method.
#'
#' @param x Numeric value(s) to evaluate the spline in.
#' @param t Numeric vector of measurement times.
#' @param y Numeric vector of values to be fitted.
#' @param gfun Choice of Greens function (choose from G1-G10).
#' @param sigma Variance, per default set to srqt(2). Effect translated into lambda.
#' @param lambda Smoothness parameter.
#' @param niter Number of iterations.
#' @param stopiter Stopping criterion.
#' @param mp Numeric value. Number of points to evaluate spline in.
#' @param constr Numeric value in (-1, 0, 1, 2).
#'               Choice of constraint.
#'               -1: No constraint.
#'               0:  Positivity contraint.
#'               1:  Monotonicity constraint.
#'               2:  Convexity/concavity constraint.
#' @param sign Numeric value in (-1, 1).
#'             1:  Increasing, convexity.
#'             -1: Decreasing, concavity.
#' @examples
#' t = generate.t(m = 25); y = generate.y()(t) + generate.noise(t, type = 2)
#' plot(t, y, col = "blue")
#' lines(t, L1spline(t, t, y, G9, lambda = 0.0001))
L1spline = function(x, t, y, gfun, sigma=sqrt(2), lambda=0.0001, niter=200, stopiter=1e-10,
                    constr=-Inf, sign=1, mp=100) {

  m    = length(t)
  beta = 1.1

  theta = psi = matrix(y, m, niter)

  L1opt = function(x) {
    out = y
    out[x > y + sqrt(2)/(2*sigma*beta)] = (-sqrt(2)/(2*sigma*beta) + x)[x > y + sqrt(2)/(2*sigma*beta)]
    out[x < y - sqrt(2)/(2*sigma*beta)] = (sqrt(2)/(2*sigma*beta) + x)[x < y - sqrt(2)/(2*sigma*beta)]
    return(out)
  }

  dist = function(x, y)
    return(t(x-y) %*% (x-y))

  if (niter > 1) {
    for (k in 1:(niter-1)) {

      beta     = beta * 1.1

      if (constr < -10)
        theta[, k+1] = thetafun(t, psi[, k], gfun, lambda=lambda/beta)
      if (constr > -10)
        theta[, k+1] = constrSpline(t, psi[, k], gfun, mp=mp, constr=constr, sign=sign, lambda=lambda)

      psi[, k+1]   = L1opt(theta[, k+1])

      if (k > 1 & dist(theta[, k+1],  theta[, k]) < stopiter) break
    }

    yout = theta[, k+1]

  } else {
    yout = y
    k    = 1
  }

  if (constr < -10)
    out = thetafun(t, yout, gfun, x = x, lambda=lambda)
  if (constr > -10)
    out = constrSpline(t, yout, gfun, mp=mp, constr=constr, sign=sign, lambda=lambda)

  attr(out, "niter") = k
  attr(out, "sigma") = sigma

  return(out)
}



