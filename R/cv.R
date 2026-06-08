#' Cross-Validate Twin-SVM Hyperparameters
#'
#' Runs k-fold cross-validation over `c1`, `c2`, and optionally `gamma`.
#'
#' @param x Numeric matrix or data frame of predictors.
#' @param y Two-class response.
#' @param c1_grid,c2_grid Positive numeric vectors.
#' @param gamma_grid Optional positive numeric vector. If `NULL`, `gamma` is
#'   left at the `tsvm()` default.
#' @param k Number of folds.
#' @param ... Additional arguments passed to [tsvm()].
#'
#' @return A `cv_tsvm` object with `best_params` and `results`.
#' @importFrom stats predict
#' @export
#'
#' @examples
#' set.seed(10)
#' dat <- gen_moons(40, noise = 0.1)
#' cv <- cv_tsvm(dat$x, dat$y, c1_grid = c(0.1, 1), c2_grid = c(0.1, 1), k = 3)
#' cv$best_params
cv_tsvm <- function(x, y, c1_grid, c2_grid, gamma_grid = NULL, k = 5, ...) {
  x <- as_numeric_matrix(x)
  y <- check_two_class_factor(y, nrow(x))
  c1_grid <- check_grid(c1_grid, "c1_grid")
  c2_grid <- check_grid(c2_grid, "c2_grid")
  if (!is.null(gamma_grid)) {
    gamma_grid <- check_grid(gamma_grid, "gamma_grid")
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2L || k > nrow(x)) {
    stop("`k` must be between 2 and `nrow(x)`.", call. = FALSE)
  }
  k <- as.integer(k)
  folds <- make_folds(y, k)
  grid <- expand.grid(
    c1 = c1_grid,
    c2 = c2_grid,
    gamma = if (is.null(gamma_grid)) NA_real_ else gamma_grid,
    KEEP.OUT.ATTRS = FALSE
  )
  extra <- list(...)
  accuracy <- numeric(nrow(grid))

  for (i in seq_len(nrow(grid))) {
    fold_acc <- numeric(k)
    for (fold in seq_len(k)) {
      test_idx <- folds[[fold]]
      train_idx <- setdiff(seq_len(nrow(x)), test_idx)
      fit_args <- c(
        list(
          x = x[train_idx, , drop = FALSE],
          y = y[train_idx],
          c1 = grid$c1[i],
          c2 = grid$c2[i]
        ),
        extra
      )
      if (!is.na(grid$gamma[i])) {
        fit_args$gamma <- grid$gamma[i]
      }
      fit <- do.call(tsvm, fit_args)
      pred <- predict(fit, x[test_idx, , drop = FALSE])
      fold_acc[fold] <- mean(pred == y[test_idx])
    }
    accuracy[i] <- mean(fold_acc)
  }
  grid$accuracy <- accuracy
  best_i <- which.max(grid$accuracy)
  best <- as.list(grid[best_i, c("c1", "c2", "gamma"), drop = FALSE])
  if (is.na(best$gamma)) {
    best$gamma <- NULL
  }

  structure(
    list(
      best_params = best,
      results = grid,
      k = k,
      call = match.call()
    ),
    class = "cv_tsvm"
  )
}

check_grid <- function(x, name) {
  if (!is.numeric(x) || length(x) < 1L || anyNA(x) || any(x <= 0)) {
    stop("`", name, "` must be a non-empty positive numeric vector.", call. = FALSE)
  }
  unique(x)
}

make_folds <- function(y, k) {
  idx <- seq_along(y)
  split(sample(idx), rep(seq_len(k), length.out = length(idx)))
}

#' Plot Twin-SVM Cross-Validation Results
#'
#' @param x A `cv_tsvm` object.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' set.seed(11)
#' dat <- gen_moons(40, noise = 0.1)
#' cv <- cv_tsvm(dat$x, dat$y, c1_grid = c(0.1, 1), c2_grid = c(0.1, 1), k = 3)
#' plot(cv)
plot.cv_tsvm <- function(x, ...) {
  res <- x$results
  has_gamma <- any(!is.na(res$gamma))
  if (length(unique(res$c1)) > 1L && length(unique(res$c2)) > 1L) {
    p <- ggplot2::ggplot(res, ggplot2::aes(.data$c1, .data$c2, fill = .data$accuracy)) +
      ggplot2::geom_tile() +
      ggplot2::geom_point(data = best_result_row(x), shape = 4, size = 3) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
      ggplot2::labs(x = "c1", y = "c2", fill = "Accuracy")
  } else {
    res$parameter <- if (length(unique(res$c1)) > 1L) res$c1 else res$c2
    p <- ggplot2::ggplot(res, ggplot2::aes(.data$parameter, .data$accuracy)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::scale_x_log10() +
      ggplot2::labs(x = "Parameter", y = "Accuracy")
  }
  if (has_gamma && length(unique(res$gamma)) > 1L) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$gamma), labeller = ggplot2::label_both)
  }
  p + .twin_theme()
}

best_result_row <- function(x) {
  x$results[which.max(x$results$accuracy), , drop = FALSE]
}
