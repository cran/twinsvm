#' Compare Twin-SVM and Standard SVM Boundaries
#'
#' Fits least-squares twin SVM, original QP twin SVM, and the standard C-SVC
#' SVM baseline on the same two-dimensional data, then plots the three decision
#' boundaries side by side.
#'
#' @param x Numeric two-column matrix or data frame.
#' @param y Two-class response.
#' @param kernel Kernel name: `"linear"`, `"rbf"`, or `"poly"`.
#' @param gamma Kernel scale.
#' @param c1,c2 Positive twin-SVM regularization parameters.
#' @param cost Positive C-SVC cost parameter.
#'
#' @return A faceted `ggplot` object.
#' @family visualization
#' @export
#'
#' @examples
#' set.seed(30)
#' dat <- gen_moons(50, noise = 0.1)
#' compare_methods(dat$x, dat$y, gamma = 1, c1 = 0.2, c2 = 0.2, cost = 1)
compare_methods <- function(x, y, kernel = "rbf", gamma = 0.5,
                            c1 = 1, c2 = 1, cost = 1) {
  x <- as_numeric_matrix(x)
  if (ncol(x) != 2L) {
    stop("`compare_methods()` requires two predictor columns.", call. = FALSE)
  }
  y <- check_two_class_factor(y, nrow(x))
  kernel <- normalize_kernel(kernel)
  gamma <- check_positive_scalar(gamma, "gamma")
  c1 <- check_positive_scalar(c1, "c1")
  c2 <- check_positive_scalar(c2, "c2")
  cost <- check_positive_scalar(cost, "cost")

  models <- list(
    ls = tsvm(x, y, method = "ls", kernel = kernel, gamma = gamma, c1 = c1, c2 = c2),
    qp = tsvm(x, y, method = "twin", kernel = kernel, gamma = gamma, c1 = c1, c2 = c2),
    svm = svms(x, y, kernel = kernel, gamma = gamma, cost = cost)
  )
  accuracy <- vapply(models, function(model) mean(stats::predict(model, x) == y), numeric(1L))
  method_labels <- c(
    sprintf("LS-TWSVM (acc %.2f)", accuracy[["ls"]]),
    sprintf("QP-TWSVM (acc %.2f)", accuracy[["qp"]]),
    sprintf("C-SVC SVM (acc %.2f)", accuracy[["svm"]])
  )

  grid_parts <- Map(function(model, label) {
    grid <- .decision_grid(model, x, n = 150, expand = 0.08)
    grid$method <- label
    grid
  }, models, method_labels)
  grid_df <- do.call(rbind, grid_parts)
  grid_df$method <- factor(grid_df$method, levels = method_labels)
  train <- data.frame(x1 = x[, 1L], x2 = x[, 2L], class = y)

  ggplot2::ggplot(grid_df, ggplot2::aes(.data$x1, .data$x2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$class), alpha = 0.22) +
    ggplot2::geom_contour(
      ggplot2::aes(z = .data$decision),
      breaks = 0,
      color = "black",
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      data = train,
      ggplot2::aes(.data$x1, .data$x2, color = .data$class),
      inherit.aes = FALSE,
      size = 1.8
    ) +
    ggplot2::facet_wrap(ggplot2::vars(.data$method), nrow = 1L) +
    ggplot2::scale_fill_manual(values = .class_colours(levels(y)), guide = "none") +
    ggplot2::scale_color_manual(values = .class_colours(levels(y)), drop = FALSE) +
    ggplot2::labs(x = "x1", y = "x2", color = "Class") +
    .twin_theme()
}
