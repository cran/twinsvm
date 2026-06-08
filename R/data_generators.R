#' Generate Two-Moons Data
#'
#' @param n Number of observations.
#' @param noise Standard deviation of Gaussian noise.
#'
#' @return A list with numeric matrix `x` and factor `y`.
#' @importFrom stats rnorm runif
#' @export
#'
#' @examples
#' dat <- gen_moons(20, noise = 0.1)
#' str(dat)
gen_moons <- function(n, noise = 0.2) {
  if (!is.numeric(n) || length(n) != 1L || n < 2L) {
    stop("`n` must be at least 2.", call. = FALSE)
  }
  n <- as.integer(n)
  n1 <- n %/% 2L
  n2 <- n - n1
  theta1 <- runif(n1, 0, pi)
  theta2 <- runif(n2, 0, pi)
  x1 <- cbind(cos(theta1), sin(theta1))
  x2 <- cbind(1 - cos(theta2), 0.5 - sin(theta2))
  x <- rbind(x1, x2)
  if (noise > 0) {
    x <- x + matrix(rnorm(length(x), sd = noise), ncol = 2)
  }
  y <- factor(c(rep("B", n1), rep("A", n2)), levels = c("B", "A"))
  list(x = x, y = y)
}

#' Generate Concentric Circles Data
#'
#' @param n Number of observations.
#' @param noise Standard deviation of Gaussian noise.
#'
#' @return A list with numeric matrix `x` and factor `y`.
#' @export
#'
#' @examples
#' dat <- gen_circles(20, noise = 0.05)
#' str(dat)
gen_circles <- function(n, noise = 0.05) {
  if (!is.numeric(n) || length(n) != 1L || n < 2L) {
    stop("`n` must be at least 2.", call. = FALSE)
  }
  n <- as.integer(n)
  n1 <- n %/% 2L
  n2 <- n - n1
  a1 <- runif(n1, 0, 2 * pi)
  a2 <- runif(n2, 0, 2 * pi)
  r1 <- 0.45 + rnorm(n1, sd = noise)
  r2 <- 1.00 + rnorm(n2, sd = noise)
  x <- rbind(cbind(r1 * cos(a1), r1 * sin(a1)),
             cbind(r2 * cos(a2), r2 * sin(a2)))
  y <- factor(c(rep("B", n1), rep("A", n2)), levels = c("B", "A"))
  list(x = x, y = y)
}
