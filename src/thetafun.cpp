// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Spline function.
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
NumericVector thetafun(NumericVector t, NumericVector y, int gfun, double lambda){
  int nrow = t.size();
  NumericMatrix G(nrow, nrow);
  for(int i = 0; i < nrow; ++i) {
    for (int j = 0; j < i+1; j++) {
      if (gfun == 1) {
        G(i, j) = G(j, i) = 1/double(6) * (- pow(t[j], 2) * pow(t[i]-1, 2) * (2*t[i]*t[j] - 3*t[i] + t[j]));
      } else if (gfun == 2) {
        G(i, j) = G(j, i) = 1/double(6) * (- t[j] * (1 - t[i]) * (pow(t[j], 2) - 2*t[i] + pow(t[i], 2)));
      } else if (gfun == 3) {
        G(i, j) = G(j, i) = 1/double(6) * (- t[j] * (3*pow(t[i], 2) + pow(t[j], 2) - 6*t[i]));
      } else if (gfun == 4) {
        G(i, j) = G(j, i) = 1/double(12) * (- t[j] * pow(t[i]-1, 2) * (t[i]*pow(t[j], 2) + 2 * pow(t[j], 2) - 3*t[i]));
      } else if (gfun == 5) {
        G(i, j) = G(j, i) = 1/double(12) * (- pow(t[j], 2) * (t[i] - 1) * (pow(t[i], 2)*t[j] - 3*pow(t[i], 2) - 2*t[i]*t[j] + 6*t[i] - 2*t[j]));
      } else if (gfun == 6) {
        G(i, j) = G(j, i) = 1/double(12) * (- pow(t[j], 2) * (3*pow(t[i], 2) - 6*t[i] + 2*t[j]));
      } else if (gfun == 7) {
        G(i, j) = G(j, i) = 1/double(12) * (pow(t[i]-1, 2) * (-3*pow(t[j], 2) + 2*t[i] + 1));
      } else if (gfun == 8) {
        G(i, j) = G(j, i) = 1/double(6) * (pow(t[i]-1, 2) * (t[i] -3*t[j] + 2));
      } else if (gfun == 9) {
        G(i, j) = G(j, i) = 1/double(6) * (pow(t[j], 2) * (3*t[i] - t[j]));
      } else if (gfun == 10) {
        G(i, j) = G(j, i) = 1/double(6) * ((t[i]-1) * (pow(t[i], 2) + 3 * pow(t[j], 2) - 2 * t[i] - 2));
      }
    }
  }
  arma::mat B;
  arma::mat C;
  arma::mat A = as<arma::mat>(G);
  arma::vec a = as<arma::vec>(y);
  C = A;
  C.diag() += lambda;
  B = arma::inv(C);
  return wrap(A * B * a);
}




