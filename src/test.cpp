// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title Spline function.
//'
//' @description Function to fit the L2 spline using some Greens function.
//'
//' @param t Numeric vector of evaluation times.
//' @param y Numeric vector to be fitted.
//' @param gfun Integer value 1-10 (choice of Greens function).
//' @param lambda Numeric (smoothness parameter).
//'
//' @return Numeric vector of fitted values.
//'
//' @author Helene Charlotte Rytgaard
//' @examples
//' t = generate.t(m = 25); y = generate.y()(t) + generate.noise(t, type = 2)
//' plot(t, y, col = "blue")
//' lines(t, thetafun(t, y, 9, lambda = 0.0001))
//' @export
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}
