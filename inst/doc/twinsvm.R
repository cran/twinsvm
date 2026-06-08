## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)

## -----------------------------------------------------------------------------
library(twinsvm)

set.seed(1)
dat <- gen_moons(100, noise = 0.12)
fit <- tsvm(dat$x, dat$y, kernel = "rbf", gamma = 2, c1 = 0.1, c2 = 0.1)
head(predict(fit, dat$x))
mean(predict(fit, dat$x) == dat$y)

## -----------------------------------------------------------------------------
plot(fit)

## -----------------------------------------------------------------------------
linear_fit <- tsvm(dat$x, dat$y, kernel = "linear")
plot(linear_fit)

## -----------------------------------------------------------------------------
cv <- cv_tsvm(
  dat$x,
  dat$y,
  c1_grid = c(0.1, 1),
  c2_grid = c(0.1, 1),
  gamma_grid = c(1, 2),
  kernel = "rbf",
  k = 3
)
cv$best_params
plot(cv)

## -----------------------------------------------------------------------------
set.seed(4)
x3 <- rbind(
  matrix(rnorm(30, -2, 0.25), ncol = 2),
  cbind(rnorm(15, 2, 0.25), rnorm(15, -2, 0.25)),
  matrix(rnorm(30, 2, 0.25), ncol = 2)
)
y3 <- factor(rep(c("alpha", "beta", "gamma"), each = 15))

multi <- tsvm(x3, y3, kernel = "linear")
head(predict(multi, x3))
head(predict(multi, x3, type = "votes"))
confusion(multi, x3, y3)

## -----------------------------------------------------------------------------
timing <- data.frame(
  n = c(40, 80, 120),
  tsvm_seconds = NA_real_,
  svms_seconds = NA_real_
)

for (i in seq_len(nrow(timing))) {
  set.seed(i)
  d <- gen_moons(timing$n[i], noise = 0.12)
  timing$tsvm_seconds[i] <- system.time(tsvm(d$x, d$y, kernel = "rbf", gamma = 2))[["elapsed"]]
  timing$svms_seconds[i] <- system.time(svms(d$x, d$y, kernel = "rbf", gamma = 2))[["elapsed"]]
}
timing

## -----------------------------------------------------------------------------
circles <- gen_circles(100, noise = 0.04)
lift_plot(circles$x, circles$y, gamma = 1)

## -----------------------------------------------------------------------------
set.seed(2)
small <- gen_moons(60, noise = 0.1)
compare_methods(small$x, small$y, gamma = 1, c1 = 0.2, c2 = 0.2, cost = 1)

## -----------------------------------------------------------------------------
anim <- morph_boundary(dat$x, dat$y, param = "gamma", range = c(0.5, 2), kernel = "rbf", n = 5)
class(anim)

