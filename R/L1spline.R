#' @title L1 spline.
#'
#' Function to compute the L1 spline, by use of the splitting method.
#'
#' @param t Numeric vector of measurement times.
#' @param y Numeric vector of values to be fitted.
#' @param y0 Numeric value to move y-values up or down. Per default y0 = 0.
#' @param gfun Choice of Greens function (1-10).
#' @param x Numeric vector. If specified, spline is computed in this set of points.
#' @param kappa Asymmetry parameter. Kappa is estimated if kappa=NULL is specified.
#'        Default is kappa = 1, corresponding to no skewness.
#' @param CI Option to compute confidence intervals.
#'                0: Nothing is computed.
#'                1: sigmahat is obtain from kappa(hat) and thetahat, and CI is computed.
#' @param lambda Smoothness parameter.
#' @param sw Integer value corresponding to number of anchor points used in sigma estimation.
#'           Default is sw=NULL, corresponding to no sigma estimation.
#' @param t0 Numeric vector of anchor points used in sigma estimation.
#'           Default is t0=NULL, so that sigma estimation is only performed when sw is not NULL.
#'           If t0 is specified, this overrules sw.
#' @param niter Number of iterations.
#' @param niterk Number of iterations for the kappa estimation.
#' @param stopiter Stopping criterion.
#' @param mp Numeric value. Number of points to evaluate spline in.
#' @param constr Integer value in (-Inf, 0, 1, 2).
#'               Choice of constraint.
#'               -Inf: No constraint.
#'               0:    Positivity contraint.
#'               1:    Monotonicity constraint.
#'               2:    Convexity/concavity constraint.
#' @param sign Numeric value in (-1, 1).
#'             1:  Increasing, convexity.
#'             -1: Decreasing, concavity.
#' @param b Numeric nonnegative value, corresponding to truncation window.
#'          If b=0, the non-adapted spline is fitted. This is default.
#' @examples
#' t = generate.t(m = 25); y = generate.y()(t) + generate.noise(t, type = 2)
#' plot(t, y, col = "blue")
#' lines(t, L1spline(t, y, 9, lambda = 0.0001))
#' @export
L1spline = function(t, y, gfun, y0=0, kappa=1, lambda=0.0001, x=Inf, niter=200,
                    sw=NULL, t0=NULL,
                    CI=0, niterk=5, stopiter=1e-10,
                    constr=-Inf, sign=1, mp=100,
                    b=0) {

  if (length(t) != length(y)) stop(paste('t and y must have same length'))
  if (!is.numeric(gfun))
    stop(paste('gfun must be in 1-10'))
  if (is.numeric(gfun) & (!(gfun %in% (1:10))))
    stop(paste('gfun must be in 1-10'))

  m = length(t)

  if (is.numeric(t0)) {
    sw = length(t0)
  }

  if (!is.numeric(kappa) | CI == 1 | is.numeric(sw)) {
    pfun = function(y, x)
      (y - x)*(y > x)
    mfun = function(y, x)
      (x - y)*(y < x)
  }

  iterfun = function(kappa, x, sigma) {

    beta = 1.1

    theta = psi = matrix(y, m, niter)

    L1optAL = function(kappa, x, sigma) {
      out = y
      if (b == 0) {
        out[x < y - kappa / (sqrt(2) * sigma * beta)] =
          (kappa / (sqrt(2) * sigma * beta) + x)[x < y - kappa / (sqrt(2) * sigma * beta)]
        out[x > 1 / (kappa * sqrt(2) * sigma * beta) + y] =
          (- 1 / (kappa * sqrt(2) * sigma * beta) +
             x)[x > 1 / (kappa * sqrt(2) * sigma * beta) + y]
      } else if (b > 0) {
        out[x - y > b]  = (y+b)[x - y > b]
        out[x - y < -b] = (y-b)[x - y < -b]
        out[x < y - kappa / (sqrt(2) * sigma * beta) - b] =
          (kappa / (sqrt(2) * sigma * beta) +
             x)[x < y - kappa / (sqrt(2) * sigma * beta) - b]
        out[x > 1 / (kappa * sqrt(2) * sigma * beta) + y + b] =
          (- 1 / (kappa * sqrt(2) * sigma * beta) +
             x)[x > 1 / (kappa * sqrt(2) * sigma * beta) + y + b]
        out[abs(x - y) < b] = x[abs(x - y) < b]
      }
      return(out)
    }

    dist = function(x, y)
      return(t(x-y) %*% (x-y))

    if (niter > 1) {
      for (k in 1:(niter-1)) {

        beta     = beta * 1.1

        if (constr < -10)
          theta[, k+1] = L2spline(t, psi[, k], gfun, y0 = y0, lambda = lambda/beta)
        if (constr > -10)
          theta[, k+1] = L2splineC(t, psi[, k], gfun, y0 = y0,
                                   mp=mp, constr=constr, sign=sign, lambda=lambda/beta)

        psi[, k+1]   = L1optAL(kappa, theta[, k+1], sigma)

        if (k > 1 & dist(theta[, k+1],  theta[, k]) < stopiter) break
      }

      yout = theta[, k+1]

    } else {

      yout = y

    }

    if (constr < -10)
      out = L2spline(t, yout, gfun, y0 = y0, lambda=lambda, x=x)
    if (constr > -10)
      out = L2splineC(t, yout, gfun, y0 = y0, mp=mp, constr=constr, sign=sign, lambda=lambda)

    attr(out, "niter") = k
    attr(out, "kappa") = kappa

    return(out)
  }

  if (!is.numeric(kappa)) {

    kappa1 = list()

    kappahat = function(theta, y, m)
      (sum(mfun(y, theta)) / sum(pfun(y, theta)))^(1/4)

    kappa1[[1]] = 1

    for (j in 1:niterk) {
      out           = iterfun(kappa1[[j]], x, 1)
      kappa1[[j+1]] = kappahat(out, y, m)
    }

    kappa = kappa1[[length(kappa1)]]
  } else {

    kappa1 = kappa

  }

  if (is.numeric(sw)) {

    if (is.numeric(t0)) {
      sigma = t0
    } else {
      qseq = seq(0, 1, length = 2*sw+1)[2:(2*sw)]
      t0   = as.numeric((quantile(t, p = qseq))[(1:length(qseq))%%2 == 1])
    }

    sigmat = function(t, t0, sigmavec) {
      if (any(sigmavec < 0)) stop('Negative variance')
      return(approx(t0, sigmavec, xout = t, rule = 2)$y)
    }

    Qs = function(theta, sigma, y, m, kappa)
      return(sum(log(sigma)) + sum(sqrt(2) * kappa / sigma * pfun(y, theta)) +
               sum(sqrt(2) / (sigma * kappa) * mfun(y, theta)))

    optfun = function(s) {
      sigma = sigmat(t, t0, abs(s))
      theta = iterfun(kappa, Inf, sigma)
      return(Qs(theta, sigma, y, m, kappa))
    }

    w0    = abs(optim(rep(1, length(t0)), optfun)$par)
    sigma = sigmat(t, t0, w0)

  } else {
    sigma = 1
  }

  out = iterfun(kappa, x, sigma)
  attr(out, "kappa")    = kappa
  attr(out, "kappavec") = unlist(kappa1)

  if (CI == 1) {

    if (!is.numeric(sw)) {
      sigmahat = function(m, y, theta, kappa) {
        return(sqrt(2) / m * (sum(pfun(y, theta)) * kappa + sum(mfun(y, theta)) / kappa))
      }
      sigmahat = sigmahat(m, y, out, kappa)
    } else {
      sigmahat = sigma
    }

    q1 = qLaplace(0.025, kappa = kappa)
    q2 = qLaplace(0.975, kappa = kappa)

    attr(out, "sigma") = sigmahat
    attr(out, "q1")    = q1
    attr(out, "q2")    = q2
    attr(out, "CI")    = sapply(c(q1, q2), function(q) out + q * sigmahat)
  }

  if (is.numeric(sw)) {
    attr(out, "sigma") = sigma
    attr(out, "t0")     = t0
    attr(out, "w0")     = w0
  }

  return(out)
}

