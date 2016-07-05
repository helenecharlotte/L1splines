################------ assuming kappa known
#' Asymmetric L1 spline.
#'
#' Function to compute the asymmetric L1 spline, by use of the continuation method.
#' Assuming either kappa (asymmetry parameter) known, or unknown (to be estimated).
#'
#' @param x Numeric value(s) to evaluate the spline in.
#' @param t Numeric vector of measurement times.
#' @param y Numeric vector of values to be fitted.
#' @param gfun Choice of Greens function (choose from G1-G10).
#' @param kappa Asymmetry parameter. Per default kappa=NULL, which means that kappa will be estimated.
#' @param lambda Smoothness parameter.
#' @param CI Option to compute confidence intervals.
#'                   0: Nothing is computed.
#'                   2: sigmahat is obtain from kappa(hat) and thetahat, and CI is computed.
#' @param niter Number of iterations.
#' @param niterk Number of iterations for the kappa estimation.
#' @param stopiter Stopping criterion.
L1splineAL = function(x, t, y, gfun, kappa=NULL, lambda=0.0001, CI = 0, niter=200, stopiter=1e-10,
                      niterk = 5) {

  m = length(t)

  iterfun = function(kappa, x) {

    beta = 1.1

    theta = psi = matrix(y, m, niter)

    L1optAL = function(kappa, x) {
      out = y
      out[x < y - kappa / (2 * beta)]     = (kappa / (2 * beta) + x)[x < y - kappa / (2 * beta)]
      out[x > 1 / (kappa * 2 * beta) + y] = (- 1 / (kappa * 2 * beta) + x)[x > 1 / (kappa * 2 * beta) + y]
      return(out)
    }

    dist = function(x, y)
      return(t(x-y) %*% (x-y))

    if (niter > 1) {
      for (k in 1:(niter-1)) {

        beta     = beta * 1.1

        theta[, k+1] = thetafun(t, psi[, k], gfun, lambda = lambda / beta)
        psi[, k+1]   = L1optAL(kappa, theta[, k+1])

        if (k > 1 & dist(theta[, k+1],  theta[, k]) < stopiter) break
      }

      yout = theta[, k+1]

    } else {

      yout = y

    }

    out = thetafun(t, yout, gfun, x = x, lambda=lambda)

    attr(out, "niter") = k
    attr(out, "kappa") = kappa

    return(out)
  }

  if (!is.numeric(kappa)) {

    kappa1 = list()

    kappahat = function(theta, y, m)
      (((sum(mfun(y, theta)))) / ((sum(pfun(y, theta)))))^(1/4)

    kappa1[[1]] = 1

    for (j in 1:niterk) {
      out           = iterfun(kappa1[[j]], x)
      kappa1[[j+1]] = kappahat(out, y, m)
    }

    kappa = kappa1[[length(kappa1)]]
  } else {

    kappa1 = kappa

  }

  out = iterfun(kappa, x)
  attr(out, "kappa")    = kappa
  attr(out, "kappavec") = unlist(kappa1)

  if (CI == 1) {
    sigmahat = function(m, y, theta, kappa) {
      pfun = function(y, x)
        (y - x)*(y > x)

      mfun = function(y, x)
        (x - y)*(y < x)

      return(sqrt(2) / m * (sum(pfun(y, theta)) * kappa + sum(mfun(y, theta)) / kappa))
    }
    sigmahat = sigmahat(m, y, out, kappa)

    q1 = qLaplace(0.025, kappa = kappa)
    q2 = qLaplace(0.975, kappa = kappa)

    attr(out, "sigma") = sigmahat
    attr(out, "q1")    = q1
    attr(out, "q2")    = q2
    attr(out, "CI")    = sapply(c(q1, q2), function(q) out + q * sigmahat)
  }

  return(out)
}
