#' Plot a Standard SVM Decision Boundary
#'
#' @param x A fitted `svms` object.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' dat <- gen_moons(40)
#' fit <- svms(dat$x, dat$y, kernel = "linear")
#' plot(fit)
plot.svms <- function(x, ...) {
  grid_df <- .decision_grid(x, x$x, n = 120, expand = 0.08)
  train <- .training_df(x)
  support <- data.frame(
    x1 = x$support_vectors[, 1L],
    x2 = x$support_vectors[, 2L]
  )
  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = grid_df,
      ggplot2::aes(.data$x1, .data$x2, fill = .data$class),
      alpha = 0.22
    ) +
    ggplot2::geom_contour(
      data = grid_df,
      ggplot2::aes(.data$x1, .data$x2, z = .data$decision),
      breaks = 0,
      color = "black",
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = train,
      ggplot2::aes(.data$x1, .data$x2, color = .data$class),
      size = 2
    ) +
    ggplot2::geom_point(
      data = support,
      ggplot2::aes(.data$x1, .data$x2),
      shape = 21,
      fill = NA,
      color = "black",
      size = 3,
      stroke = 0.8
    ) +
    ggplot2::scale_fill_manual(values = .class_colours(x$levels), guide = "none") +
    ggplot2::scale_color_manual(values = .class_colours(x$levels), drop = FALSE) +
    ggplot2::labs(x = "x1", y = "x2", color = "Class") +
    .twin_theme()
}

#' Plot a Twin-SVM Decision Boundary
#'
#' @param x A fitted `tsvm` object.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @export
#'
#' @examples
#' dat <- gen_moons(40)
#' fit <- tsvm(dat$x, dat$y)
#' plot(fit)
plot.tsvm <- function(x, ...) {
  grid_df <- .decision_grid(x, x$x, n = 120, expand = 0.08)
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = grid_df,
      ggplot2::aes(.data$x1, .data$x2, fill = .data$class),
      alpha = 0.22
    ) +
    ggplot2::geom_contour(
      data = grid_df,
      ggplot2::aes(.data$x1, .data$x2, z = .data$decision),
      breaks = 0,
      color = "black",
      linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = .training_df(x),
      ggplot2::aes(.data$x1, .data$x2, color = .data$class),
      size = 2
    ) +
    ggplot2::scale_fill_manual(values = .class_colours(x$levels), guide = "none") +
    ggplot2::scale_color_manual(values = .class_colours(x$levels), drop = FALSE) +
    ggplot2::labs(x = "x1", y = "x2", color = "Class") +
    .twin_theme()

  if (x$kernel == "linear") {
    p <- add_twin_planes(p, x)
  }
  p
}

add_twin_planes <- function(p, x) {
  if (abs(x$w1[2L]) > 1e-12) {
    p <- p + ggplot2::geom_abline(
      slope = -x$w1[1L] / x$w1[2L],
      intercept = -x$b1 / x$w1[2L],
      color = "#dc2626",
      linetype = "dashed"
    )
  }
  if (abs(x$w2[2L]) > 1e-12) {
    p <- p + ggplot2::geom_abline(
      slope = -x$w2[1L] / x$w2[2L],
      intercept = -x$b2 / x$w2[2L],
      color = "#2563eb",
      linetype = "dashed"
    )
  }
  p
}
