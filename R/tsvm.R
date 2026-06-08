#' Fit a Twin Support Vector Machine
#'
#' Fits a twin support vector machine. With two classes this uses the validated
#' binary twin-SVM path: level 1 of `y` is class B, level 2 is class A, plane 1
#' is close to class A, and plane 2 is close to class B. With three or more
#' classes, the function fits one binary twin SVM for each class pair and
#' predicts by majority vote. Multiclass ties are resolved by choosing the class
#' that appears first in the factor level order.
#'
#' @param x Numeric matrix or data frame of predictors.
#' @param y Response with at least two classes.
#' @param method Twin-SVM method. `"ls"` fits least-squares twin SVM;
#'   `"twin"` fits the original box-constrained dual formulation.
#' @param kernel Kernel name.
#' @param c1,c2 Positive regularization parameters.
#' @param gamma Kernel scale. Defaults to `1 / ncol(x)`.
#' @param degree Polynomial degree.
#' @param coef0 Polynomial offset.
#' @param eps Ridge term added to every linear solve.
#'
#' @return A fitted `tsvm` object for two classes, or a `tsvm_multiclass` object
#'   for three or more classes.
#' @export
#'
#' @references
#' Jayadeva, Khemchandani, R., and Chandra, S. (2007). Twin support vector
#' machines for pattern classification. *IEEE Transactions on Pattern Analysis
#' and Machine Intelligence*, 29(5), 905-910.
#'
#' Kumar, M. A. and Gopal, M. (2009). Least squares twin support vector
#' machines for pattern classification. *Expert Systems with Applications*,
#' 36(4), 7535-7543.
#'
#' @examples
#' set.seed(3)
#' dat <- gen_moons(50, noise = 0.05)
#' fit <- tsvm(dat$x, dat$y)
#' predict(fit, dat$x[1:4, ])
tsvm <- function(x, y, method = c("ls", "twin"),
                 kernel = c("linear", "rbf", "poly"), c1 = 1, c2 = 1,
                 gamma = NULL, degree = 3, coef0 = 1, eps = 1e-6) {
  x <- as_numeric_matrix(x)
  y <- check_class_factor(y, nrow(x))
  if (nlevels(y) > 2L) {
    method <- match.arg(method)
    kernel <- normalize_kernel(kernel)
    c1 <- check_positive_scalar(c1, "c1")
    c2 <- check_positive_scalar(c2, "c2")
    gamma <- check_positive_scalar(default_gamma(gamma, x), "gamma")
    degree <- check_positive_scalar(degree, "degree")
    coef0 <- check_nonnegative_scalar(coef0, "coef0")
    eps <- check_positive_scalar(eps, "eps")
    binary_fitter <- function(x_pair, y_pair) {
      tsvm(
        x_pair, y_pair,
        method = method,
        kernel = kernel,
        c1 = c1,
        c2 = c2,
        gamma = gamma,
        degree = degree,
        coef0 = coef0,
        eps = eps
      )
    }
    return(.ovo_fit(
      binary_fitter = binary_fitter,
      x = x,
      y = y,
      object_class = "tsvm_multiclass",
      call = match.call(),
      extra = list(
        method = method,
        kernel = kernel,
        c1 = c1,
        c2 = c2,
        gamma = gamma,
        degree = as.integer(degree),
        coef0 = coef0,
        eps = eps
      )
    ))
  }
  y <- check_two_class_factor(y, nrow(x))
  method <- match.arg(method)
  kernel <- normalize_kernel(kernel)
  c1 <- check_positive_scalar(c1, "c1")
  c2 <- check_positive_scalar(c2, "c2")
  gamma <- check_positive_scalar(default_gamma(gamma, x), "gamma")
  degree <- check_positive_scalar(degree, "degree")
  coef0 <- check_nonnegative_scalar(coef0, "coef0")
  eps <- check_positive_scalar(eps, "eps")
  parts <- class_split(x, y)
  if (kernel == "linear") {
    if (method == "ls") {
      fit <- lstsvm_linear_cpp(parts$A, parts$B, c1, c2, eps)
    } else {
      fit <- qptsvm_linear_cpp(parts$A, parts$B, c1, c2, eps)
    }
    model_fields <- list(
      w1 = fit$w1,
      b1 = fit$b1,
      w2 = fit$w2,
      b2 = fit$b2,
      norm1 = fit$norm1,
      norm2 = fit$norm2
    )
  } else {
    basis <- rbind(parts$A, parts$B)
    if (method == "ls") {
      fit <- lstsvm_kernel_cpp(
        parts$A, parts$B, basis, kernel, gamma, as.integer(degree), coef0,
        c1, c2, eps
      )
    } else {
      fit <- qptsvm_kernel_cpp(
        parts$A, parts$B, basis, kernel, gamma, as.integer(degree), coef0,
        c1, c2, eps
      )
    }
    model_fields <- list(
      basis = basis,
      u1 = fit$u1,
      b1 = fit$b1,
      u2 = fit$u2,
      b2 = fit$b2,
      norm1 = fit$norm1,
      norm2 = fit$norm2
    )
  }

  object <- c(
    list(
      x = x,
      y = y,
      levels = levels(y),
      method = method,
      kernel = kernel,
      c1 = c1,
      c2 = c2,
      gamma = gamma,
      degree = as.integer(degree),
      coef0 = coef0,
      eps = eps,
      n_features = ncol(x),
      call = match.call()
    ),
    model_fields
  )
  structure(object, class = "tsvm")
}

#' @export
print.tsvm <- function(x, ...) {
  cat("Twin support vector machine\n")
  cat("  Method:", x$method, "\n")
  cat("  Kernel:", x$kernel, "\n")
  cat("  Classes: B =", x$levels[1L], ", A =", x$levels[2L], "\n")
  invisible(x)
}

#' Extract Twin-SVM Coefficients
#'
#' @param object A fitted `tsvm` object.
#' @param ... Unused.
#'
#' @return A list with the two plane coefficients.
#' @export
#'
#' @examples
#' dat <- gen_moons(30)
#' fit <- tsvm(dat$x, dat$y)
#' coef(fit)
coef.tsvm <- function(object, ...) {
  if (object$kernel == "linear") {
    list(w1 = object$w1, b1 = object$b1, w2 = object$w2, b2 = object$b2)
  } else {
    list(u1 = object$u1, b1 = object$b1, u2 = object$u2, b2 = object$b2)
  }
}

#' @export
predict.tsvm <- function(object, newdata, decision.values = FALSE, ...) {
  x <- check_newdata(object, newdata)
  if (object$kernel == "linear") {
    f1 <- drop(x %*% object$w1 + object$b1)
    f2 <- drop(x %*% object$w2 + object$b2)
  } else {
    k <- kernel_matrix_cpp(
      x, object$basis, object$kernel, object$gamma,
      object$degree, object$coef0
    )
    f1 <- drop(k %*% object$u1 + object$b1)
    f2 <- drop(k %*% object$u2 + object$b2)
  }
  d1 <- abs(f1) / object$norm1
  d2 <- abs(f2) / object$norm2
  decision <- d2 - d1
  if (decision.values) {
    return(decision)
  }
  factor(ifelse(d1 <= d2, object$levels[2L], object$levels[1L]), levels = object$levels)
}
