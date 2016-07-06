##--- specifying green's functions

#' Greens function 1.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G1 = function(s, t)
  1/6 * (t-s)^3 * (t > s) - 1/6 * t^2 * (s-1)^2 * (2*s*t - 3*s + t)

#' Greens function 2.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G2 = function(s, t)
  1/6 * (t-s)^3 * (t > s) - 1/6 * t * (1-s) * (t^2 - 2*s + s^2)

#' Greens function 3.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G3 = function(s, t)
  1/6 * (t-s)^3 * (t > s) - 1/6 * t * (3*s^2 + t^2 - 6*s)

#' Greens function 4.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G4 = function(s, t)
  1/6 * (t-s)^3 * (t > s) - 1/12 * t*(s-1)^2*(s*t^2 + 2*t^2-3*s)

#' Greens function 5.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G5 = function(s, t)
  1/6 * (t-s)^3 * (t > s) - 1/12 * t^2 * (s-1) * (s^2*t - 3*s^2 - 2*s*t + 6*s - 2*t)

#' Greens function 6.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G6 = function(s, t)
  1/6 * (t-s)^3 * (t > s) - 1/12 * t^2 * (3*s^2 - 6*s + 2*t)

#' Greens function 7.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G7 = function(s, t)
  1/6 * (t-s)^3 * (t > s) + 1/12 * (s-1)^2 * (-3*t^2 + 2*s + 1)

#' Greens function 8.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G8 = function(s, t)
  1/6 * (t-s)^3 * (t > s) + 1/6 * (s-1)^2 * (s - 3*t + 2)

#' Greens function 9.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G9 = function(s, t)
  1/6 * (t-s)^3 * (t > s) + 1/6*t^2 * (3*s - t)

#' Greens function 10.
#'
#' Function to compute the Greens function in two time points.
#'
#' @param s Numeric value.
#' @param t Numeric value.
G10 = function(s, t)
  1/6 * (t-s)^3 * (t > s) + 1/6 * (s-1) * (s^2 + 3*t^2 - 2*s - 2)

##--- theta fit function

#' Spline function.
#'
#' Function to fit the L2 spline using some Greens function.
#'
#' @param t Numeric vector of observation times.
#' @param y Numeric vector to be fitted.
#' @param gfun Choice of Greens function (1-10).
#' @param x Numeric vector. If specified, spline is computed in this set of points.
#' @param lambda Smoothness parameter.
#' @examples
#' t = generate.t(m = 25); y = generate.y()(t) + generate.noise(t, type = 2)
#' plot(t, y, col = "blue")
#' lines(t, thetafun(t, y, G9, lambda = 0.0001))
#' @export
L2spline = function(t, y, gfun, lambda=0.001, x = Inf) {
  m = length(t)

  if (length(t) != length(y)) stop(paste('t and y must have same length'))
  if (!is.numeric(gfun))
    stop(paste('gfun must be in 1-10'))
  if (is.numeric(gfun) & (!(gfun %in% (1:10))))
    stop(paste('gfun must be in 1-10'))

  if (any(x < Inf))
    out = outerG(t, x, gfun, 0) %*% thetafunG(t, y, gfun, lambda)
  else
    out = thetafun(t, y, gfun, lambda)

  return(as.numeric(out))
}

