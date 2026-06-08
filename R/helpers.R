as_numeric_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  } else {
    x <- as.matrix(x)
  }
  storage.mode(x) <- "double"
  if (!is.numeric(x)) {
    stop("`", name, "` must be numeric.", call. = FALSE)
  }
  if (anyNA(x)) {
    stop("`", name, "` must not contain missing values.", call. = FALSE)
  }
  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop("`", name, "` must have at least one row and one column.", call. = FALSE)
  }
  x
}

check_two_class_factor <- function(y, n = NULL) {
  if (!is.factor(y)) {
    y <- factor(y)
  }
  if (!is.null(n) && length(y) != n) {
    stop("`y` must have the same length as `nrow(x)`.", call. = FALSE)
  }
  y <- droplevels(y)
  if (nlevels(y) != 2L) {
    stop("`y` must have exactly two classes.", call. = FALSE)
  }
  y
}

check_class_factor <- function(y, n = NULL) {
  if (!is.factor(y)) {
    y <- factor(y)
  }
  if (!is.null(n) && length(y) != n) {
    stop("`y` must have the same length as `nrow(x)`.", call. = FALSE)
  }
  y <- droplevels(y)
  if (nlevels(y) < 2L) {
    stop("`y` must have at least two classes.", call. = FALSE)
  }
  y
}

normalize_kernel <- function(kernel) {
  match.arg(kernel, c("linear", "rbf", "poly"))
}

default_gamma <- function(gamma, x) {
  if (is.null(gamma)) {
    1 / ncol(x)
  } else {
    gamma
  }
}

check_positive_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value <= 0) {
    stop("`", name, "` must be a positive numeric scalar.", call. = FALSE)
  }
  value
}

check_nonnegative_scalar <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) || value < 0) {
    stop("`", name, "` must be a non-negative numeric scalar.", call. = FALSE)
  }
  value
}

label_from_decision <- function(decision, levels) {
  factor(ifelse(decision >= 0, levels[2L], levels[1L]), levels = levels)
}

check_newdata <- function(object, newdata) {
  x <- as_numeric_matrix(newdata, "newdata")
  if (ncol(x) != object$n_features) {
    stop("`newdata` must have ", object$n_features, " columns.", call. = FALSE)
  }
  x
}

class_split <- function(x, y) {
  list(
    B = x[y == levels(y)[1L], , drop = FALSE],
    A = x[y == levels(y)[2L], , drop = FALSE]
  )
}
