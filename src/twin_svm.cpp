// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static arma::mat augment_with_ones(const arma::mat& X) {
  return arma::join_horiz(X, arma::ones<arma::vec>(X.n_rows));
}

static arma::mat kernel_matrix_internal(const arma::mat& X, const arma::mat& Z,
                                        const std::string& kernel,
                                        double gamma, int degree,
                                        double coef0) {
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

// [[Rcpp::export]]
Rcpp::List lstsvm_linear_cpp(const arma::mat& A, const arma::mat& B,
                             double c1, double c2, double eps) {
  arma::mat H = augment_with_ones(A);
  arma::mat G = augment_with_ones(B);
  arma::uword p = H.n_cols;
  arma::mat I = arma::eye<arma::mat>(p, p);

  arma::mat M1 = H.t() * H + c1 * (G.t() * G) + eps * I;
  arma::vec r1 = -c1 * G.t() * arma::ones<arma::vec>(G.n_rows);
  arma::vec z1 = arma::solve(M1, r1, arma::solve_opts::likely_sympd);

  arma::mat M2 = G.t() * G + c2 * (H.t() * H) + eps * I;
  arma::vec r2 = c2 * H.t() * arma::ones<arma::vec>(H.n_rows);
  arma::vec z2 = arma::solve(M2, r2, arma::solve_opts::likely_sympd);

  arma::vec w1 = z1.head(p - 1);
  arma::vec w2 = z2.head(p - 1);
  double norm1 = std::sqrt(std::max(arma::dot(w1, w1), 1e-30));
  double norm2 = std::sqrt(std::max(arma::dot(w2, w2), 1e-30));

  return Rcpp::List::create(
    Rcpp::Named("w1") = w1,
    Rcpp::Named("b1") = z1(p - 1),
    Rcpp::Named("w2") = w2,
    Rcpp::Named("b2") = z2(p - 1),
    Rcpp::Named("norm1") = norm1,
    Rcpp::Named("norm2") = norm2
  );
}

// [[Rcpp::export]]
Rcpp::List lstsvm_kernel_cpp(const arma::mat& A, const arma::mat& B,
                             const arma::mat& C, const std::string& kernel,
                             double gamma, int degree, double coef0,
                             double c1, double c2, double eps) {
  arma::mat H = augment_with_ones(kernel_matrix_internal(A, C, kernel, gamma, degree, coef0));
  arma::mat G = augment_with_ones(kernel_matrix_internal(B, C, kernel, gamma, degree, coef0));
  arma::uword p = H.n_cols;
  arma::mat I = arma::eye<arma::mat>(p, p);

  arma::mat M1 = H.t() * H + c1 * (G.t() * G) + eps * I;
  arma::vec r1 = -c1 * G.t() * arma::ones<arma::vec>(G.n_rows);
  arma::vec z1 = arma::solve(M1, r1, arma::solve_opts::likely_sympd);

  arma::mat M2 = G.t() * G + c2 * (H.t() * H) + eps * I;
  arma::vec r2 = c2 * H.t() * arma::ones<arma::vec>(H.n_rows);
  arma::vec z2 = arma::solve(M2, r2, arma::solve_opts::likely_sympd);

  arma::vec u1 = z1.head(p - 1);
  arma::vec u2 = z2.head(p - 1);
  arma::mat KC = kernel_matrix_internal(C, C, kernel, gamma, degree, coef0);
  double norm1 = std::sqrt(std::max(arma::as_scalar(u1.t() * KC * u1), 1e-30));
  double norm2 = std::sqrt(std::max(arma::as_scalar(u2.t() * KC * u2), 1e-30));

  return Rcpp::List::create(
    Rcpp::Named("u1") = u1,
    Rcpp::Named("b1") = z1(p - 1),
    Rcpp::Named("u2") = u2,
    Rcpp::Named("b2") = z2(p - 1),
    Rcpp::Named("norm1") = norm1,
    Rcpp::Named("norm2") = norm2
  );
}
