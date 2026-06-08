.twin_class_colours <- c(B = "#F8766D", A = "#00B3B3")

.class_colours <- function(class_levels) {
  if (length(class_levels) != 2L) {
    stop("Visualization helpers require exactly two class levels.", call. = FALSE)
  }
  structure(unname(.twin_class_colours), names = class_levels)
}

.twin_theme <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

.decision_grid <- function(model, x, n = 200, expand = 0.1) {
  x <- as_numeric_matrix(x)
  if (ncol(x) != 2L) {
    stop("Decision-boundary plots require exactly two predictor columns.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 2L) {
    stop("`n` must be at least 2.", call. = FALSE)
  }
  expand <- check_nonnegative_scalar(expand, "expand")
  n <- as.integer(n)

  rx <- range(x[, 1L])
  ry <- range(x[, 2L])
  padx <- diff(rx) * expand
  pady <- diff(ry) * expand
  xs <- seq(rx[1L] - padx, rx[2L] + padx, length.out = n)
  ys <- seq(ry[1L] - pady, ry[2L] + pady, length.out = n)
  grid <- as.matrix(expand.grid(x1 = xs, x2 = ys, KEEP.OUT.ATTRS = FALSE))
  decision <- stats::predict(model, grid, decision.values = TRUE)

  data.frame(
    x1 = grid[, 1L],
    x2 = grid[, 2L],
    class = label_from_decision(decision, model$levels),
    decision = decision
  )
}

.training_df <- function(object) {
  data.frame(
    x1 = object$x[, 1L],
    x2 = object$x[, 2L],
    class = object$y
  )
}
