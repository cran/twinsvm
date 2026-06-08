#' Static RBF Kernel Lift Plot
#'
#' Maps two-dimensional data to a third RBF height and draws a static oblique
#' projection. The height is `exp(-gamma * ||x - center||^2)`, so points near
#' the center rise while points farther away stay low. A translucent plane is
#' drawn halfway between the two class mean heights.
#'
#' @param x Numeric two-column matrix or data frame.
#' @param y Two-class response.
#' @param gamma Positive RBF scale.
#' @param center Optional numeric length-two center for the RBF bump. If `NULL`,
#'   the centroid of the inner class is used.
#'
#' @return A `ggplot` object.
#' @family visualization
#' @export
#'
#' @examples
#' set.seed(20)
#' dat <- gen_circles(80, noise = 0.04)
#' lift_plot(dat$x, dat$y, gamma = 1)
lift_plot <- function(x, y, gamma = 1, center = NULL) {
  lift <- .lift_geometry(x, y, gamma, center, caller = "lift_plot")
  points <- .project_lift(lift$data$x1, lift$data$x2, lift$data$z)
  base <- .project_lift(lift$data$x1, lift$data$x2, 0)
  plot_data <- data.frame(
    lift$data,
    sx = points$sx,
    sy = points$sy,
    sx0 = base$sx,
    sy0 = base$sy
  )
  plane <- .lift_plane_polygon(lift)

  ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = plane,
      ggplot2::aes(.data$sx, .data$sy),
      fill = "grey55",
      color = "grey35",
      alpha = 0.22,
      linewidth = 0.35
    ) +
    ggplot2::geom_segment(
      data = plot_data,
      ggplot2::aes(
        x = .data$sx0,
        y = .data$sy0,
        xend = .data$sx,
        yend = .data$sy
      ),
      color = "grey65",
      alpha = 0.35,
      linewidth = 0.25
    ) +
    ggplot2::geom_point(
      data = plot_data,
      ggplot2::aes(.data$sx, .data$sy, color = .data$class),
      size = 2.2,
      alpha = 0.9
    ) +
    ggplot2::coord_equal() +
    ggplot2::scale_color_manual(values = .class_colours(levels(lift$data$class)), drop = FALSE) +
    ggplot2::labs(
      x = "projected x1/x2",
      y = "projected RBF height",
      color = "Class"
    ) +
    .twin_theme()
}

#' Interactive RBF Kernel Lift Plot
#'
#' Draws the same RBF lift as [lift_plot()] as a rotatable `plotly` 3D chart.
#' This function is optional; the package does not require `plotly` to build or
#' check.
#'
#' @param x Numeric two-column matrix or data frame.
#' @param y Two-class response.
#' @param gamma Positive RBF scale.
#' @param center Optional numeric length-two center for the RBF bump. If `NULL`,
#'   the centroid of the inner class is used.
#'
#' @return A `plotly` object.
#' @family visualization
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(21)
#' dat <- gen_circles(80, noise = 0.04)
#' if (requireNamespace("plotly", quietly = TRUE)) {
#'   lift_plotly(dat$x, dat$y, gamma = 1)
#' }
#' }
lift_plotly <- function(x, y, gamma = 1, center = NULL) {
  if (!.plotly_available()) {
    stop(
      "`lift_plotly()` requires the optional plotly package. ",
      "Install it with install.packages(\"plotly\"), or use `lift_plot()` ",
      "for the static ggplot version.",
      call. = FALSE
    )
  }
  lift <- .lift_geometry(x, y, gamma, center, caller = "lift_plotly")
  x_seq <- seq(lift$x_range[1L], lift$x_range[2L], length.out = 2L)
  y_seq <- seq(lift$y_range[1L], lift$y_range[2L], length.out = 2L)
  plane_z <- matrix(lift$plane_z, nrow = 2L, ncol = 2L)

  p <- plotly::plot_ly()
  p <- plotly::add_trace(
    p,
    data = lift$data,
    x = ~x1,
    y = ~x2,
    z = ~z,
    color = ~class,
    colors = .class_colours(levels(lift$data$class)),
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 4),
    name = "points"
  )
  p <- plotly::add_surface(
    p,
    x = x_seq,
    y = y_seq,
    z = plane_z,
    opacity = 0.25,
    showscale = FALSE,
    name = "separating plane"
  )
  plotly::layout(
    p,
    scene = list(
      xaxis = list(title = "x1"),
      yaxis = list(title = "x2"),
      zaxis = list(title = "RBF height")
    )
  )
}

#' Compatibility Alias for the RBF Kernel Lift Plot
#'
#' `kernel_lift()` is kept for existing code. New code should call
#' [lift_plot()].
#'
#' @param x Numeric two-column matrix or data frame.
#' @param y Two-class response.
#' @param gamma Positive RBF scale.
#' @param center Optional numeric length-two center for the RBF bump. If `NULL`,
#'   the centroid of the inner class is used.
#'
#' @return A `ggplot` object.
#' @family visualization
#' @export
#'
#' @examples
#' set.seed(22)
#' dat <- gen_circles(60, noise = 0.04)
#' kernel_lift(dat$x, dat$y, gamma = 1)
kernel_lift <- function(x, y, gamma = 1, center = NULL) {
  message("`kernel_lift()` is kept for compatibility; use `lift_plot()` for new code.")
  lift_plot(x, y, gamma = gamma, center = center)
}

.lift_geometry <- function(x, y, gamma, center, caller) {
  x <- as_numeric_matrix(x)
  if (ncol(x) != 2L) {
    stop("`", caller, "()` requires two predictor columns.", call. = FALSE)
  }
  y <- check_two_class_factor(y, nrow(x))
  gamma <- check_positive_scalar(gamma, "gamma")
  if (is.null(center)) {
    overall_center <- colMeans(x)
    centered <- sweep(x, 2L, overall_center, "-")
    class_radius <- tapply(rowSums(centered^2), y, mean)
    inner_level <- names(which.min(class_radius))
    center <- colMeans(x[y == inner_level, , drop = FALSE])
  } else {
    if (!is.numeric(center) || length(center) != 2L || anyNA(center)) {
      stop("`center` must be a numeric vector of length two.", call. = FALSE)
    }
    center <- as.numeric(center)
  }
  centered <- sweep(x, 2L, center, "-")
  z <- exp(-gamma * rowSums(centered^2))
  mean_heights <- tapply(z, y, mean)

  list(
    data = data.frame(x1 = x[, 1L], x2 = x[, 2L], z = z, class = y),
    center = center,
    plane_z = mean(mean_heights),
    x_range = range(x[, 1L]),
    y_range = range(x[, 2L])
  )
}

.project_lift <- function(x1, x2, z) {
  data.frame(
    sx = x1 + 0.55 * x2 * cos(pi / 6),
    sy = 1.4 * z + 0.55 * x2 * sin(pi / 6)
  )
}

.lift_plane_polygon <- function(lift) {
  corners <- data.frame(
    x1 = c(lift$x_range[1L], lift$x_range[2L], lift$x_range[2L], lift$x_range[1L]),
    x2 = c(lift$y_range[1L], lift$y_range[1L], lift$y_range[2L], lift$y_range[2L]),
    z = lift$plane_z
  )
  .project_lift(corners$x1, corners$x2, corners$z)
}

.plotly_available <- function() {
  !isTRUE(getOption("twinsvm.force_no_plotly")) &&
    requireNamespace("plotly", quietly = TRUE)
}
