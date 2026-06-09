# twinsvm

`twinsvm` brings twin support vector machines to R with a standard C-SVC SVM baseline, one-vs-one multiclass classification, and visualization-first helpers for small to moderate problems.

## Install locally

```r
install.packages("twinsvm_0.0.2.tar.gz", repos = NULL, type = "source")
library(twinsvm)
```

## Quick start

```r
set.seed(1)
dat <- gen_moons(100, noise = 0.12)

fit <- tsvm(dat$x, dat$y, kernel = "rbf", gamma = 2)
pred <- predict(fit, dat$x)
mean(pred == dat$y)

plot(fit)
```

## Multiclass

`tsvm()` and `svms()` support three or more classes with one-vs-one majority voting. Ties are resolved by the first factor level.

```r
set.seed(4)
x3 <- rbind(
  matrix(rnorm(30, -2, 0.25), ncol = 2),
  cbind(rnorm(15, 2, 0.25), rnorm(15, -2, 0.25)),
  matrix(rnorm(30, 2, 0.25), ncol = 2)
)
y3 <- factor(rep(c("alpha", "beta", "gamma"), each = 15))

multi <- tsvm(x3, y3, kernel = "linear")
predict(multi, x3[1:5, , drop = FALSE])
predict(multi, x3[1:5, , drop = FALSE], type = "votes")
confusion(multi, x3, y3)
```

## Visualization

Show the RBF lift as a static pseudo-3D ggplot:

```r
circles <- gen_circles(100, noise = 0.04)
lift_plot(circles$x, circles$y, gamma = 1)
```

Compare least-squares twin SVM, original QP twin SVM, and the standard SVM baseline side by side:

```r
compare_methods(dat$x, dat$y, gamma = 1, c1 = 0.2, c2 = 0.2, cost = 1)
```

Compare against the standard SVM baseline:

```r
base <- svms(dat$x, dat$y, kernel = "rbf", cost = 1, gamma = 2)
mean(predict(base, dat$x) == dat$y)
```

Tune the twin SVM:

```r
cv <- cv_tsvm(
  dat$x, dat$y,
  c1_grid = c(0.1, 1),
  c2_grid = c(0.1, 1),
  gamma_grid = c(1, 2),
  kernel = "rbf",
  k = 3
)
cv$best_params
plot(cv)
```

## Validation note

For the standard SVM baseline, tests compare predictions against `e1071`, the LIBSVM-backed R package. For twin SVM there is no existing R implementation to compare against, so tests check the defining plane-distance behavior, nonlinear kernel improvement, and agreement between least-squares and original QP twin-SVM variants. Kernel twin-SVM models invert an `(n + 1)` matrix and are intended for small to moderate data sets.
