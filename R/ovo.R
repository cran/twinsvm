.ovo_fit <- function(binary_fitter, x, y, object_class, call, extra = list()) {
  class_levels <- levels(y)
  pair_indices <- utils::combn(seq_along(class_levels), 2L)
  n_pairs <- ncol(pair_indices)
  models <- vector("list", n_pairs)
  pair_names <- vector("list", n_pairs)

  for (i in seq_len(n_pairs)) {
    pair <- pair_indices[, i]
    names_i <- class_levels[pair]
    rows <- y %in% names_i
    y_pair <- factor(y[rows], levels = names_i)
    counts <- tabulate(as.integer(y_pair), nbins = 2L)
    if (any(counts < 2L)) {
      stop(
        "Cannot fit one-vs-one pair `", names_i[1L], "` vs `", names_i[2L],
        "`: each class must have at least two observations.",
        call. = FALSE
      )
    }
    models[[i]] <- tryCatch(
      binary_fitter(x[rows, , drop = FALSE], y_pair),
      error = function(e) {
        stop(
          "Failed to fit one-vs-one pair `", names_i[1L], "` vs `", names_i[2L],
          "`: ", conditionMessage(e),
          call. = FALSE
        )
      }
    )
    pair_names[[i]] <- names_i
  }

  pair_df <- data.frame(
    class1 = vapply(pair_names, `[`, character(1L), 1L),
    class2 = vapply(pair_names, `[`, character(1L), 2L),
    stringsAsFactors = FALSE
  )

  object <- c(
    list(
      models = models,
      pair_indices = t(pair_indices),
      pairs = pair_df,
      levels = class_levels,
      k = length(class_levels),
      n_models = n_pairs,
      x = x,
      y = y,
      n_features = ncol(x),
      call = call
    ),
    extra
  )
  structure(object, class = object_class)
}

.ovo_predict <- function(object, newdata, type = c("class", "votes"), ...) {
  dots <- list(...)
  if ("decision.values" %in% names(dots)) {
    stop(
      "`decision.values` is only defined for binary models; use `type = \"votes\"` ",
      "for multiclass diagnostics.",
      call. = FALSE
    )
  }
  type <- match.arg(type)
  x <- check_newdata(object, newdata)
  votes <- matrix(
    0L,
    nrow = nrow(x),
    ncol = object$k,
    dimnames = list(NULL, object$levels)
  )

  rows <- seq_len(nrow(x))
  for (i in seq_along(object$models)) {
    pred <- stats::predict(object$models[[i]], x)
    cols <- match(as.character(pred), object$levels)
    votes[cbind(rows, cols)] <- votes[cbind(rows, cols)] + 1L
  }

  if (type == "votes") {
    return(votes)
  }
  factor(object$levels[max.col(votes, ties.method = "first")], levels = object$levels)
}

#' Predict from a Multiclass Twin SVM
#'
#' Predicts from a one-vs-one multiclass twin SVM. Each binary model votes for
#' one class. Ties are resolved deterministically by choosing the class that
#' appears first in the training factor level order.
#'
#' @param object A fitted `tsvm_multiclass` object.
#' @param newdata Numeric matrix or data frame.
#' @param type Output type. `"class"` returns predicted class labels;
#'   `"votes"` returns the vote matrix.
#' @param ... Unused. `decision.values` is not supported for multiclass
#'   objects because OVO decision values do not have a single unambiguous scale.
#'
#' @return A factor of predicted classes, or an integer vote matrix when
#'   `type = "votes"`.
#' @family multiclass
#' @export
#'
#' @examples
#' set.seed(40)
#' x <- rbind(
#'   matrix(rnorm(20, -2, 0.2), ncol = 2),
#'   matrix(rnorm(20, 0, 0.2), ncol = 2),
#'   matrix(rnorm(20, 2, 0.2), ncol = 2)
#' )
#' y <- factor(rep(c("a", "b", "c"), each = 10))
#' fit <- tsvm(x, y, kernel = "linear")
#' predict(fit, x[1:3, , drop = FALSE])
predict.tsvm_multiclass <- function(object, newdata, type = c("class", "votes"), ...) {
  .ovo_predict(object, newdata, type = type, ...)
}

#' Print a Multiclass Twin SVM
#'
#' @param x A fitted `tsvm_multiclass` object.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @family multiclass
#' @export
#'
#' @examples
#' set.seed(41)
#' x <- rbind(
#'   matrix(rnorm(20, -2, 0.2), ncol = 2),
#'   matrix(rnorm(20, 0, 0.2), ncol = 2),
#'   matrix(rnorm(20, 2, 0.2), ncol = 2)
#' )
#' y <- factor(rep(c("a", "b", "c"), each = 10))
#' print(tsvm(x, y, kernel = "linear"))
print.tsvm_multiclass <- function(x, ...) {
  cat("One-vs-one twin support vector machine\n")
  cat("  Method:", x$method, "\n")
  cat("  Kernel:", x$kernel, "\n")
  cat("  Classes:", x$k, "(", paste(x$levels, collapse = " / "), ")\n", sep = "")
  cat("  Binary models:", x$n_models, "\n")
  invisible(x)
}

#' Predict from a Multiclass Standard SVM
#'
#' Predicts from a one-vs-one multiclass standard SVM. Each binary model votes
#' for one class. Ties are resolved deterministically by choosing the class that
#' appears first in the training factor level order.
#'
#' @param object A fitted `svms_multiclass` object.
#' @param newdata Numeric matrix or data frame.
#' @param type Output type. `"class"` returns predicted class labels;
#'   `"votes"` returns the vote matrix.
#' @param ... Unused. `decision.values` is not supported for multiclass
#'   objects because OVO decision values do not have a single unambiguous scale.
#'
#' @return A factor of predicted classes, or an integer vote matrix when
#'   `type = "votes"`.
#' @family multiclass
#' @export
#'
#' @examples
#' set.seed(45)
#' x <- rbind(
#'   matrix(rnorm(8, -2, 0.2), ncol = 2),
#'   matrix(rnorm(8, 0, 0.2), ncol = 2),
#'   matrix(rnorm(8, 2, 0.2), ncol = 2)
#' )
#' y <- factor(rep(c("a", "b", "c"), each = 4))
#' fit <- svms(x, y, kernel = "linear", max_passes = 2, max_iter = 100)
#' predict(fit, x[1:3, , drop = FALSE])
predict.svms_multiclass <- function(object, newdata, type = c("class", "votes"), ...) {
  .ovo_predict(object, newdata, type = type, ...)
}

#' Print a Multiclass Standard SVM
#'
#' @param x A fitted `svms_multiclass` object.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#' @family multiclass
#' @export
#'
#' @examples
#' set.seed(46)
#' x <- rbind(
#'   matrix(rnorm(8, -2, 0.2), ncol = 2),
#'   matrix(rnorm(8, 0, 0.2), ncol = 2),
#'   matrix(rnorm(8, 2, 0.2), ncol = 2)
#' )
#' y <- factor(rep(c("a", "b", "c"), each = 4))
#' print(svms(x, y, kernel = "linear", max_passes = 2, max_iter = 100))
print.svms_multiclass <- function(x, ...) {
  cat("One-vs-one C-SVC support vector machine\n")
  cat("  Kernel:", x$kernel, "\n")
  cat("  Classes:", x$k, "(", paste(x$levels, collapse = " / "), ")\n", sep = "")
  cat("  Binary models:", x$n_models, "\n")
  invisible(x)
}
