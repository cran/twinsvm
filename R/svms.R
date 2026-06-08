#' Fit a Standard C-SVC Support Vector Machine
#'
#' Fits a C-SVC support vector machine with a Platt SMO solver. With two
#' classes this uses the validated binary path. With three or more classes, the
#' function fits one binary SVM for each class pair and predicts by majority
#' vote. Multiclass ties are resolved by choosing the class that appears first
#' in the factor level order.
#'
#' @param x Numeric matrix or data frame of predictors.
#' @param y Response with at least two classes. In the binary path, level 1 is
#'   the negative class and level 2 is the positive class.
#' @param kernel Kernel name: `"linear"`, `"rbf"`, or `"poly"`.
#' @param cost Positive C-SVC cost parameter.
#' @param gamma Kernel scale. Defaults to `1 / ncol(x)`.
#' @param degree Polynomial degree.
#' @param coef0 Polynomial offset.
#' @param tol SMO tolerance.
#' @param max_passes Maximum consecutive passes without alpha changes.
#' @param max_iter Maximum SMO iterations.
#'
#' @return A fitted `svms` object for two classes, or a `svms_multiclass` object
#'   for three or more classes.
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- gen_moons(40, noise = 0.1)
#' fit <- svms(dat$x, dat$y, kernel = "linear", cost = 1)
#' predict(fit, dat$x[1:3, ])
svms <- function(x, y, kernel = c("linear", "rbf", "poly"), cost = 1,
                 gamma = NULL, degree = 3, coef0 = 1, tol = 1e-3,
                 max_passes = 10L, max_iter = 10000L) {
  x <- as_numeric_matrix(x)
  y <- check_class_factor(y, nrow(x))
  if (nlevels(y) > 2L) {
    kernel <- normalize_kernel(kernel)
    cost <- check_positive_scalar(cost, "cost")
    gamma <- check_positive_scalar(default_gamma(gamma, x), "gamma")
    degree <- check_positive_scalar(degree, "degree")
    coef0 <- check_nonnegative_scalar(coef0, "coef0")
    binary_fitter <- function(x_pair, y_pair) {
      svms(
        x_pair, y_pair,
        kernel = kernel,
        cost = cost,
        gamma = gamma,
        degree = degree,
        coef0 = coef0,
        tol = tol,
        max_passes = max_passes,
        max_iter = max_iter
      )
    }
    return(.ovo_fit(
      binary_fitter = binary_fitter,
      x = x,
      y = y,
      object_class = "svms_multiclass",
      call = match.call(),
      extra = list(
        kernel = kernel,
        cost = cost,
        gamma = gamma,
        degree = as.integer(degree),
        coef0 = coef0,
        tol = tol,
        max_passes = as.integer(max_passes),
        max_iter = as.integer(max_iter)
      )
    ))
  }
  y <- check_two_class_factor(y, nrow(x))
  kernel <- normalize_kernel(kernel)
  cost <- check_positive_scalar(cost, "cost")
  gamma <- check_positive_scalar(default_gamma(gamma, x), "gamma")
  degree <- check_positive_scalar(degree, "degree")
  coef0 <- check_nonnegative_scalar(coef0, "coef0")
  y_num <- ifelse(y == levels(y)[2L], 1, -1)

  fit <- smo_cpp(
    x, y_num, cost, kernel, gamma, as.integer(degree), coef0,
    tol, as.integer(max_passes), as.integer(max_iter)
  )

  structure(
    list(
      x = x,
      y = y,
      levels = levels(y),
      kernel = kernel,
      cost = cost,
      gamma = gamma,
      degree = as.integer(degree),
      coef0 = coef0,
      alpha = fit$alpha,
      b = fit$b,
      support_indices = fit$support_indices,
      support_vectors = fit$support_vectors,
      support_alpha = fit$support_alpha,
      support_y = fit$support_y,
      n_support = length(fit$support_indices),
      n_features = ncol(x),
      call = match.call()
    ),
    class = "svms"
  )
}

#' @export
print.svms <- function(x, ...) {
  cat("C-SVC support vector machine\n")
  cat("  Kernel:", x$kernel, "\n")
  cat("  Classes:", paste(x$levels, collapse = " / "), "\n")
  cat("  Support vectors:", x$n_support, "\n")
  invisible(x)
}
