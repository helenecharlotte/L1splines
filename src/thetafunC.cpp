// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' @title Outer function for gfun.
//'
//' @description Function to compute the matrix of Green's function evaluation for any two vectors t and tp.
//'
//' @param t Numeric vector of observation times.
//' @param tp Numeric vector evaluation times.
//' @param gfun Integer value 1-10 (choice of Greens function).
//' @param inv Logical (want inverse matrix or not)
//'
//' @return Numeric vector of fitted values.
//'
//' @author Helene Charlotte Rytgaard
//' @export
// [[Rcpp::export]]
NumericMatrix outerG(NumericVector t, NumericVector tp, int gfun, bool inv){
  int ncol = t.size();
  int nrow = tp.size();

  NumericMatrix G1(nrow, ncol);
  for(int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; j++) {
      if (gfun == 1) {
        G1(i, j) = 1/double(6) * (- pow(t[j], 2) * pow(tp[i]-1, 2) * (2*tp[i]*t[j] - 3*tp[i] + t[j]));
      } else if (gfun == 2) {
        G1(i, j) = 1/double(6) * (- t[j] * (1 - tp[i]) * (pow(t[j], 2) - 2*tp[i] + pow(tp[i], 2)));
      } else if (gfun == 3) {
        G1(i, j) = 1/double(6) * (- t[j] * (3*pow(tp[i], 2) + pow(t[j], 2) - 6*tp[i]));
      } else if (gfun == 4) {
        G1(i, j) =  1/double(12) * (- t[j] * pow(tp[i]-1, 2) * (tp[i]*pow(t[j], 2) + 2 * pow(t[j], 2) - 3*tp[i]));
      } else if (gfun == 5) {
        G1(i, j) = 1/double(12) * (- pow(t[j], 2) * (tp[i] - 1) * (pow(tp[i], 2)*t[j] - 3*pow(tp[i], 2) - 2*tp[i]*t[j] + 6*tp[i] - 2*t[j]));
      } else if (gfun == 6) {
        G1(i, j) = 1/double(12) * (- pow(t[j], 2) * (3*pow(tp[i], 2) - 6*tp[i] + 2*t[j]));
      } else if (gfun == 7) {
        G1(i, j) = 1/double(12) * (pow(tp[i]-1, 2) * (-3*pow(t[j], 2) + 2*tp[i] + 1));
      } else if (gfun == 8) {
        G1(i, j) = 1/double(6) * (pow(tp[i]-1, 2) * (tp[i] -3*t[j] + 2));
      } else if (gfun == 9) {
        G1(i, j) = 1/double(6) * (pow(t[j], 2) * (3*tp[i] - t[j]));
      } else if (gfun == 10) {
        G1(i, j) = 1/double(6) * ((tp[i]-1) * (pow(tp[i], 2) + 3 * pow(t[j], 2) - 2 * tp[i] - 2));
      }
      if (t[j] > tp[i]) {
        G1(i, j) += 1/double(6) * pow(t[j]-tp[i], 3);
      }
    }
  }

  if (inv) {
    arma::mat G  = as<arma::mat>(G1);
    arma::mat G2 = arma::pinv(G);
    return wrap(G2);
  } else {
    return wrap(G1);
  }
}

//' @title General solution, Green's spline.
//'
//' @description Function that takes the observations and computes the general spline solution. Must be combined with outerG() to obtain solution in some set of points.
//'
//' @param t Numeric vector of observation times.
//' @param y Numeric vector of observations.
//' @param gfun Integer value 1-10 (choice of Greens function).
//' @param lambda Numeric (smoothness parameter).
//'
//' @return Numeric vector of fitted values.
//'
//' @author Helene Charlotte Rytgaard
//' @export
// [[Rcpp::export]]
NumericMatrix thetafunG(NumericVector t, NumericVector y, int gfun, double lambda){
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
  return wrap(B * a);
}











NumericVector thetafunC(NumericVector t, NumericVector tp, NumericVector y, int gfun, double lambda){
  int nrow = t.size();
  int mp   = tp.size();

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

  NumericMatrix G1(nrow, mp);
  for(int i = 0; i < nrow; ++i) {
    for (int j = 0; j < mp; j++) {
      if (gfun == 1) {
        G1(i, j) = 1/double(6) * (- pow(tp[j], 2) * pow(t[i]-1, 2) * (2*t[i]*tp[j] - 3*t[i] + tp[j]));
      } else if (gfun == 2) {
        G1(i, j) = 1/double(6) * (- tp[j] * (1 - t[i]) * (pow(tp[j], 2) - 2*t[i] + pow(t[i], 2)));
      } else if (gfun == 3) {
        G1(i, j) = 1/double(6) * (- tp[j] * (3*pow(t[i], 2) + pow(tp[j], 2) - 6*t[i]));
      } else if (gfun == 4) {
        G1(i, j) =  1/double(12) * (- tp[j] * pow(t[i]-1, 2) * (t[i]*pow(tp[j], 2) + 2 * pow(tp[j], 2) - 3*t[i]));
      } else if (gfun == 5) {
        G1(i, j) = 1/double(12) * (- pow(tp[j], 2) * (t[i] - 1) * (pow(t[i], 2)*tp[j] - 3*pow(t[i], 2) - 2*t[i]*tp[j] + 6*t[i] - 2*tp[j]));
      } else if (gfun == 6) {
        G1(i, j) = 1/double(12) * (- pow(tp[j], 2) * (3*pow(t[i], 2) - 6*t[i] + 2*tp[j]));
      } else if (gfun == 7) {
        G1(i, j) = 1/double(12) * (pow(t[i]-1, 2) * (-3*pow(tp[j], 2) + 2*t[i] + 1));
      } else if (gfun == 8) {
        G1(i, j) = 1/double(6) * (pow(t[i]-1, 2) * (t[i] -3*tp[j] + 2));
      } else if (gfun == 9) {
        G1(i, j) = 1/double(6) * (pow(tp[j], 2) * (3*t[i] - tp[j]));
      } else if (gfun == 10) {
        G1(i, j) = 1/double(6) * ((t[i]-1) * (pow(t[i], 2) + 3 * pow(tp[j], 2) - 2 * t[i] - 2));
      }
      if (tp[j] > t[i]) {
        G1(i, j) += 1/double(6) * pow(tp[j]-t[i], 3);
      }
    }
  }

  NumericMatrix Gp(mp, mp);
  for(int i = 0; i < mp; ++i) {
    for (int j = 0; j < i+1; j++) {
      if (gfun == 1) {
        Gp(i, j) = Gp(j, i) = 1/double(6) * (- pow(t[j], 2) * pow(t[i]-1, 2) * (2*t[i]*t[j] - 3*t[i] + t[j]));
      } else if (gfun == 2) {
        Gp(i, j) = Gp(j, i) = 1/double(6) * (- t[j] * (1 - t[i]) * (pow(t[j], 2) - 2*t[i] + pow(t[i], 2)));
      } else if (gfun == 3) {
        Gp(i, j) = Gp(j, i) = 1/double(6) * (- t[j] * (3*pow(t[i], 2) + pow(t[j], 2) - 6*t[i]));
      } else if (gfun == 4) {
        Gp(i, j) = Gp(j, i) = 1/double(12) * (- t[j] * pow(t[i]-1, 2) * (t[i]*pow(t[j], 2) + 2 * pow(t[j], 2) - 3*t[i]));
      } else if (gfun == 5) {
        Gp(i, j) = Gp(j, i) = 1/double(12) * (- pow(t[j], 2) * (t[i] - 1) * (pow(t[i], 2)*t[j] - 3*pow(t[i], 2) - 2*t[i]*t[j] + 6*t[i] - 2*t[j]));
      } else if (gfun == 6) {
        Gp(i, j) = Gp(j, i) = 1/double(12) * (- pow(t[j], 2) * (3*pow(t[i], 2) - 6*t[i] + 2*t[j]));
      } else if (gfun == 7) {
        Gp(i, j) = Gp(j, i) = 1/double(12) * (pow(t[i]-1, 2) * (-3*pow(t[j], 2) + 2*t[i] + 1));
      } else if (gfun == 8) {
        Gp(i, j) = Gp(j, i) = 1/double(6) * (pow(t[i]-1, 2) * (t[i] -3*t[j] + 2));
      } else if (gfun == 9) {
        Gp(i, j) = Gp(j, i) = 1/double(6) * (pow(t[j], 2) * (3*t[i] - t[j]));
      } else if (gfun == 10) {
        Gp(i, j) = Gp(j, i) = 1/double(6) * ((t[i]-1) * (pow(t[i], 2) + 3 * pow(t[j], 2) - 2 * t[i] - 2));
      }
    }
  }

  arma::mat D;
  arma::vec eta;

  arma::mat B;
  arma::mat C;
  arma::mat A = as<arma::mat>(G);
  arma::vec a = as<arma::vec>(y);
  C = A;
  C.diag() += lambda;
  B = arma::inv(C);

  arma::mat A1 = as<arma::mat>(G1);
  arma::mat Ap = as<arma::mat>(Gp);
  arma::mat Bp = arma::pinv(Ap);

  D = lambda * Bp;
  D.diag() += 1;
  eta = A1 * B * a;

  //List out;
  //out["Dmat"] = wrap(D);
  //out["eta"]  = wrap(eta);

  return wrap(D);
}




