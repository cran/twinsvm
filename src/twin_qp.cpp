// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static arma::mat augment_with_ones_qp(const arma::mat& X) {
  return arma::join_horiz(X, arma::ones<arma::vec>(X.n_rows));
}

static arma::mat kernel_matrix_qp(const arma::mat& X, const arma::mat& Z,
                                  const std::string& kernel,
                                  double gamma, int degree, double coef0) {
  if (kernel == "linear") {
    return X * Z.t();
  }
  if (kernel == "poly") {
    return arma::pow(gamma * (X * Z.t()) + coef0, degree);
  }
  if (kernel == "rbf") {
    arma::vec x2 = arma::sum(arma::square(X), 1);
    arma::vec z2 = arma::sum(arma::square(Z), 1);
    arma::mat dist2 = arma::repmat(x2, 1, Z.n_rows) +
      arma::repmat(z2.t(), X.n_rows, 1) - 2.0 * (X * Z.t());
    return arma::exp(-gamma * dist2);
  }
  Rcpp::stop("Unknown kernel.");
}

static arma::vec box_coordinate_ascent(const arma::mat& Q, double c,
                                       int max_iter, double tol) {
  arma::uword n = Q.n_rows;
  arma::vec alpha(n, arma::fill::zeros);
  arma::vec Qa(n, arma::fill::zeros);

  for (int iter = 0; iter < max_iter; ++iter) {
    double max_delta = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      double qii = std::max(Q(i, i), 1e-12);
      double old_value = alpha(i);
      double new_value = old_value + (1.0 - Qa(i)) / qii;
      if (new_value < 0.0) new_value = 0.0;
      if (new_value > c) new_value = c;
      double delta = new_value - old_value;
      if (std::abs(delta) > 0.0) {
        alpha(i) = new_value;
        Qa += delta * Q.col(i);
        max_delta = std::max(max_delta, std::abs(delta));
      }
    }
    if (max_delta < tol) {
      break;
    }
  }
  return alpha;
}

static Rcpp::List qp_from_designs(const arma::mat& H, const arma::mat& G,
                                  double c1, double c2, double eps,
                                  int max_iter, double tol) {
  arma::uword p = H.n_cols;
  arma::mat I = arma::eye<arma::mat>(p, p);
  arma::mat HH = H.t() * H + eps * I;
  arma::mat GG = G.t() * G + eps * I;

  arma::mat invHH_Gt = arma::solve(HH, G.t(), arma::solve_opts::likely_sympd);
  arma::mat Q1 = G * invHH_Gt;
  Q1 = 0.5 * (Q1 + Q1.t());
  arma::vec a = box_coordinate_ascent(Q1, c1, max_iter, tol);
  arma::vec z1 = -arma::solve(HH, G.t() * a, arma::solve_opts::likely_sympd);

  arma::mat invGG_Ht = arma::solve(GG, H.t(), arma::solve_opts::likely_sympd);
  arma::mat Q2 = H * invGG_Ht;
  Q2 = 0.5 * (Q2 + Q2.t());
  arma::vec b = box_coordinate_ascent(Q2, c2, max_iter, tol);
  arma::vec z2 = arma::solve(GG, H.t() * b, arma::solve_opts::likely_sympd);

  return Rcpp::List::create(
    Rcpp::Named("z1") = z1,
    Rcpp::Named("z2") = z2,
    Rcpp::Named("dual1") = a,
    Rcpp::Named("dual2") = b
  );
}

// [[Rcpp::export]]
Rcpp::List qptsvm_linear_cpp(const arma::mat& A, const arma::mat& B,
                             double c1, double c2, double eps,
                             int max_iter = 10000, double tol = 1e-6) {
  arma::mat H = augment_with_ones_qp(A);
  arma::mat G = augment_with_ones_qp(B);
  Rcpp::List raw = qp_from_designs(H, G, c1, c2, eps, max_iter, tol);
  arma::vec z1 = raw["z1"];
  arma::vec z2 = raw["z2"];
  arma::uword p = z1.n_elem;
  arma::vec w1 = z1.head(p - 1);
  arma::vec w2 = z2.head(p - 1);

  return Rcpp::List::create(
    Rcpp::Named("w1") = w1,
    Rcpp::Named("b1") = z1(p - 1),
    Rcpp::Named("w2") = w2,
    Rcpp::Named("b2") = z2(p - 1),
    Rcpp::Named("norm1") = std::sqrt(std::max(arma::dot(w1, w1), 1e-30)),
    Rcpp::Named("norm2") = std::sqrt(std::max(arma::dot(w2, w2), 1e-30)),
    Rcpp::Named("dual1") = raw["dual1"],
    Rcpp::Named("dual2") = raw["dual2"]
  );
}

// [[Rcpp::export]]
Rcpp::List qptsvm_kernel_cpp(const arma::mat& A, const arma::mat& B,
                             const arma::mat& C, const std::string& kernel,
                             double gamma, int degree, double coef0,
                             double c1, double c2, double eps,
                             int max_iter = 10000, double tol = 1e-6) {
  arma::mat H = augment_with_ones_qp(kernel_matrix_qp(A, C, kernel, gamma, degree, coef0));
  arma::mat G = augment_with_ones_qp(kernel_matrix_qp(B, C, kernel, gamma, degree, coef0));
  Rcpp::List raw = qp_from_designs(H, G, c1, c2, eps, max_iter, tol);
  arma::vec z1 = raw["z1"];
  arma::vec z2 = raw["z2"];
  arma::uword p = z1.n_elem;
  arma::vec u1 = z1.head(p - 1);
  arma::vec u2 = z2.head(p - 1);
  arma::mat KC = kernel_matrix_qp(C, C, kernel, gamma, degree, coef0);

  return Rcpp::List::create(
    Rcpp::Named("u1") = u1,
    Rcpp::Named("b1") = z1(p - 1),
    Rcpp::Named("u2") = u2,
    Rcpp::Named("b2") = z2(p - 1),
    Rcpp::Named("norm1") = std::sqrt(std::max(arma::as_scalar(u1.t() * KC * u1), 1e-30)),
    Rcpp::Named("norm2") = std::sqrt(std::max(arma::as_scalar(u2.t() * KC * u2), 1e-30)),
    Rcpp::Named("dual1") = raw["dual1"],
    Rcpp::Named("dual2") = raw["dual2"]
  );
}
