// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
double plusfun(double x) {
  if (x > 0) {
    return x;
  } else {
    return 0;
  }
}

// [[Rcpp::export]]
double G1(double s, double t) {
  double out =  1/(double)6 * (pow(t-s, 3) * plusfun(t-s) - pow(t, 2) * pow(s-1, 2) * (2*s*t-3*s+t));
  return out;
}


//' something something
//'
//' @description
//' @param x A matrix.
//' @return A matrix of same size as x.
//' @author Helene Charlotte Rytgaard
//' @examples
//' @export
// [[Rcpp::export]]
NumericMatrix outerfun(NumericVector t, NumericVector y, int gfun, double lambda){
  int nrow = t.size();
  NumericMatrix out(nrow, nrow);
  for(int i = 0; i < nrow; ++i) {
    out(i, i) = 1/double(6) * (- pow(t[i], 2) * pow(t[i]-1, 2) * (2*t[i]*t[i] - 3*t[i] + t[i])) + lambda;
    for (int j = 0; j < i; j++) {
      if (gfun == 1) {
        out(i, j) = out(j, i) = 1/double(6) * (- pow(t[j], 2) * pow(t[i]-1, 2) * (2*t[i]*t[j] - 3*t[i] + t[j]));
      }
    }
  }
  arma::mat A = as<arma::mat>(out);
  arma::mat B = pinv(A);
  return wrap(B);
}




