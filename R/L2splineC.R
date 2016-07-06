#' @title Constrained L2 spline.
#'
#' Function to compute a constrained L2 spline.
#'
#' @param t Numeric vector of measurement times.
#' @param y Numeric vector of values to be fitted.
#' @param gfun Choice of Greens function (choose from G1-G10).
#' @param mp Integer value. Number of points to evaluate spline in.
#' @param constr Integer value in (0, 1, 2).
#'               Choice of constraint.
#'               0:  Positivity contraint.
#'               1:  Monotonicity constraint.
#'               2:  Convexity/concavity constraint.
#' @param sign Numeric value in (-1, 1).
#'             1:  Increasing, convexity.
#'             -1: Decreasing, concavity.
#' @param lambda Smoothness parameter.
#' @examples
#' t = generate.t(m = 25); y = generate.y(type = 2)(t) + generate.noise(t) + generate.noise(t)
#' plot(t, y, col = "blue")
#' S = L2splineC(t, y, 3, constr=1)
#' lines(attr(S, "t"), S)
#' t = generate.t(m = 25); y = -generate.y(type = 3)(t) + generate.noise(t) + generate.noise(t)
#' plot(t, y, col = "blue")
#' S = L2splineC(t, y, 3, constr=2, sign=-1)
#' lines(attr(S, "t"), S)
#' t = generate.t(m = 40); y = generate.y(type = 4)(t) + generate.noise(t)
#' plot(t, y, col = "blue")
#' S = L2splineC(t, y, 1, constr=0)
#' lines(t, S)
#' @export
L2splineC = function(t, y, gfun, mp=100, constr=0, sign=1, lambda=0.0001) {

  if (length(t) != length(y)) stop(paste('t and y must have same length'))
  if (!is.numeric(gfun))
    stop(paste('gfun must be in 1-10'))
  if (is.numeric(gfun) & (!(gfun %in% (1:10))))
    stop(paste('gfun must be in 1-10'))
  if (!is.numeric(constr))
    stop(paste('constr must be in 0-2'))
  if (is.numeric(constr) & (!(constr %in% (0:2))))
    stop(paste('constr must be in 0-2'))
  if (!is.numeric(sign))
    stop(paste('sign must be in (-1, 1)'))
  if (is.numeric(sign) & (!(sign %in% c(-1, 1))))
    stop(paste('sign must be in (-1, 1)'))

  tp = seq(0, 1, length = mp)

  thetahat  = outerG(t, tp, gfun, 0) %*% thetafunG(t, y, gfun, lambda)
  L         = lambda * outerG(tp, tp, 1, 1)

  Dmat = diag(mp) + L
  eta = Dmat %*% thetahat

  if (constr == 1) {
    Amat = diag(x = -1, mp)
    Amat[cbind(1:(mp - 1), 2:mp)] = 1
    Amat = t(Amat)
    Amat = Amat[, -mp]
  }

  if (constr == 0) {
    Amat = diag(mp)
  }

  if (constr == 2) {
    Amat = diag(x = 0.5, mp)
    Amat[cbind(1:(mp - 1), 2:mp)] = -1
    Amat[cbind(1:(mp - 2), 3:mp)] = 0.5
    Amat = t(Amat)
    Amat = Amat[, -c(mp-1, mp)]
  }

  proj       = quadprog::solve.QP(Dmat, eta, sign*Amat)
  thetahat_c = proj$solution

  out = spline(tp, thetahat_c, xout = t)$y

  return(out)
}
