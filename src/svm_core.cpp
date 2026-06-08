// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

static arma::mat kernel_matrix_internal(const arma::mat& X, const arma::mat& Z,
                                        const std::string& kernel,
                                        double gamma, int degree,
                                        double coef0) {
  arma::mat K(X.n_rows, Z.n_rows, arma::fill::zeros);
  if (kernel == "linear") {
    K = X * Z.t();
  } else if (kernel == "poly") {
    K = arma::pow(gamma * (X * Z.t()) + coef0, degree);
  } else if (kernel == "rbf") {
    arma::vec x2 = arma::sum(arma::square(X), 1);
    arma::vec z2 = arma::sum(arma::square(Z), 1);
    K = arma::repmat(x2, 1, Z.n_rows) +
      arma::repmat(z2.t(), X.n_rows, 1) - 2.0 * (X * Z.t());
    K = arma::exp(-gamma * K);
  } else {
    Rcpp::stop("Unknown kernel.");
  }
  return K;
}

// [[Rcpp::export]]
arma::mat kernel_matrix_cpp(const arma::mat& X, const arma::mat& Z,
                            const std::string& kernel, double gamma,
                            int degree, double coef0) {
  return kernel_matrix_internal(X, Z, kernel, gamma, degree, coef0);
}

static double svm_decision_one(const arma::vec& alpha, const arma::vec& y,
                               const arma::mat& K, double b, arma::uword i) {
  return arma::dot(alpha % y, K.col(i)) + b;
}

// [[Rcpp::export]]
Rcpp::List smo_cpp(const arma::mat& X, const arma::vec& y, double cost,
                   const std::string& kernel, double gamma, int degree,
                   double coef0, double tol, int max_passes, int max_iter) {
  const arma::uword n = X.n_rows;
  arma::mat K = kernel_matrix_internal(X, X, kernel, gamma, degree, coef0);
  arma::vec alpha(n, arma::fill::zeros);
  double b = 0.0;
  int passes = 0;
  int iter = 0;

  while (passes < max_passes && iter < max_iter) {
    int changed = 0;
    for (arma::uword i = 0; i < n; ++i) {
      double Ei = svm_decision_one(alpha, y, K, b, i) - y(i);
      bool violates = (y(i) * Ei < -tol && alpha(i) < cost) ||
        (y(i) * Ei > tol && alpha(i) > 0.0);
      if (!violates) {
        continue;
      }

      arma::uword j = i;
      if (n > 1) {
        j = static_cast<arma::uword>(std::floor(R::runif(0.0, static_cast<double>(n - 1))));
        if (j >= i) {
          j += 1;
        }
      }

      double Ej = svm_decision_one(alpha, y, K, b, j) - y(j);
      double ai_old = alpha(i);
      double aj_old = alpha(j);

      double L, H;
      if (y(i) != y(j)) {
        L = std::max(0.0, aj_old - ai_old);
        H = std::min(cost, cost + aj_old - ai_old);
      } else {
        L = std::max(0.0, ai_old + aj_old - cost);
        H = std::min(cost, ai_old + aj_old);
      }
      if (std::abs(L - H) < 1e-12) {
        continue;
      }

      double eta = 2.0 * K(i, j) - K(i, i) - K(j, j);
      if (eta >= 0.0) {
        continue;
      }

      alpha(j) = aj_old - y(j) * (Ei - Ej) / eta;
      if (alpha(j) > H) alpha(j) = H;
      if (alpha(j) < L) alpha(j) = L;
      if (std::abs(alpha(j) - aj_old) < 1e-5) {
        alpha(j) = aj_old;
        continue;
      }

      alpha(i) = ai_old + y(i) * y(j) * (aj_old - alpha(j));

      double b1 = b - Ei -
        y(i) * (alpha(i) - ai_old) * K(i, i) -
        y(j) * (alpha(j) - aj_old) * K(i, j);
      double b2 = b - Ej -
        y(i) * (alpha(i) - ai_old) * K(i, j) -
        y(j) * (alpha(j) - aj_old) * K(j, j);

      if (alpha(i) > 0.0 && alpha(i) < cost) {
        b = b1;
      } else if (alpha(j) > 0.0 && alpha(j) < cost) {
        b = b2;
      } else {
        b = 0.5 * (b1 + b2);
      }
      changed += 1;
    }
    passes = changed == 0 ? passes + 1 : 0;
    iter += 1;
  }

  arma::uvec sv = arma::find(alpha > 1e-7);
  arma::uvec support_indices = sv + 1;
  arma::mat support_vectors = X.rows(sv);
  arma::vec support_alpha = alpha.elem(sv);
  arma::vec support_y = y.elem(sv);
  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("b") = b,
    Rcpp::Named("support_indices") = support_indices,
    Rcpp::Named("support_vectors") = support_vectors,
    Rcpp::Named("support_alpha") = support_alpha,
    Rcpp::Named("support_y") = support_y,
    Rcpp::Named("iterations") = iter
  );
}
