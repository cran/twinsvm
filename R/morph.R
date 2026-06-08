#' Animate a Decision Boundary Across a Hyperparameter Range
#'
#' Refits a model over a sequence of hyperparameter values and returns a
#' `gganimate` object showing the boundary change.
#'
#' @param x Numeric two-column matrix or data frame.
#' @param y Two-class response.
#' @param param Hyperparameter to vary.
#' @param range Numeric length-two range for the hyperparameter.
#' @param model Model family: `"tsvm"` or `"svms"`.
#' @param n Number of frames.
#' @param ... Additional arguments passed to the model fit function.
#'
#' @return A `gganim` object.
#' @export
#'
#' @examples
#' \donttest{
#' if (interactive()) {
#'   dat <- gen_moons(40)
#'   morph_boundary(dat$x, dat$y, param = "gamma", range = c(0.5, 2), n = 4)
#' }
#' }
morph_boundary <- function(x, y, param = c("gamma", "cost", "c1"),
                           range, model = c("tsvm", "svms"), n = 30, ...) {
  if (!requireNamespace("gganimate", quietly = TRUE)) {
    stop("Package `gganimate` is required for `morph_boundary()`.", call. = FALSE)
  }
  x <- as_numeric_matrix(x)
  if (ncol(x) != 2L) {
    stop("`morph_boundary()` requires two predictor columns.", call. = FALSE)
  }
  y <- check_two_class_factor(y, nrow(x))
  param <- match.arg(param)
  model <- match.arg(model)
  if (!is.numeric(range) || length(range) != 2L || anyNA(range) || any(range <= 0)) {
    stop("`range` must contain two positive numeric values.", call. = FALSE)
  }
  if (model == "tsvm" && param == "cost") {
    stop("Use `param = \"c1\"` for `tsvm`; `cost` belongs to `svms`.", call. = FALSE)
  }
  if (model == "svms" && param == "c1") {
    stop("Use `param = \"cost\"` for `svms`; `c1` belongs to `tsvm`.", call. = FALSE)
  }
  values <- seq(range[1L], range[2L], length.out = as.integer(n))
  extra <- list(...)
  frames <- vector("list", length(values))

  for (i in seq_along(values)) {
    fit_args <- c(list(x = x, y = y), extra)
    fit_args[[param]] <- values[i]
    fit <- if (model == "tsvm") {
      do.call(tsvm, fit_args)
    } else {
      do.call(svms, fit_args)
    }
    frame <- .decision_grid(fit, x, n = 80)
    frame$value <- values[i]
    frames[[i]] <- frame
  }
  grid_df <- do.call(rbind, frames)
  train <- data.frame(x1 = x[, 1L], x2 = x[, 2L], class = y)

  ggplot2::ggplot(grid_df, ggplot2::aes(.data$x1, .data$x2)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$class), alpha = 0.22) +
    ggplot2::geom_contour(ggplot2::aes(z = .data$decision), breaks = 0, color = "black") +
    ggplot2::geom_point(
      data = train,
      ggplot2::aes(.data$x1, .data$x2, color = .data$class),
      inherit.aes = FALSE,
      size = 2
    ) +
    ggplot2::scale_fill_manual(values = .class_colours(levels(y)), guide = "none") +
    ggplot2::scale_color_manual(values = .class_colours(levels(y)), drop = FALSE) +
    ggplot2::labs(x = "x1", y = "x2", color = "Class") +
    .twin_theme() +
    gganimate::transition_states(.data$value, transition_length = 1, state_length = 1) +
    ggplot2::labs(title = paste0(param, " = {closest_state}"))
}
