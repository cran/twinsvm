test_that("standard SVM agrees with e1071 on moons", {
  skip_if_not_installed("e1071")
  set.seed(12)
  dat <- gen_moons(120, noise = 0.1)

  fit_rbf <- svms(dat$x, dat$y, kernel = "rbf", cost = 1, gamma = 1)
  ref_rbf <- e1071::svm(
    x = dat$x,
    y = dat$y,
    kernel = "radial",
    cost = 1,
    gamma = 1,
    scale = FALSE
  )
  expect_gte(mean(predict(fit_rbf, dat$x) == predict(ref_rbf, dat$x)), 0.99)
  expect_lte(abs(fit_rbf$n_support - ref_rbf$tot.nSV), 1)

  fit_linear <- svms(dat$x, dat$y, kernel = "linear", cost = 1)
  ref_linear <- e1071::svm(
    x = dat$x,
    y = dat$y,
    kernel = "linear",
    cost = 1,
    scale = FALSE
  )
  expect_gte(mean(predict(fit_linear, dat$x) == predict(ref_linear, dat$x)), 0.99)
  expect_lte(abs(fit_linear$n_support - ref_linear$tot.nSV), 2)
})

test_that("validated baseline numbers stay stable", {
  set.seed(12)
  dat <- gen_moons(120, noise = 0.1)
  fit_rbf <- svms(dat$x, dat$y, kernel = "rbf", cost = 1, gamma = 1)
  expect_equal(fit_rbf$n_support, 25)
  expect_equal(mean(predict(fit_rbf, dat$x) == dat$y), 1)

  x_b <- matrix(rnorm(80, mean = -2, sd = 0.25), ncol = 2)
  x_a <- matrix(rnorm(80, mean = 2, sd = 0.25), ncol = 2)
  x <- rbind(x_b, x_a)
  y <- factor(c(rep("B", 40), rep("A", 40)), levels = c("B", "A"))
  fit <- tsvm(x, y)
  expect_gt(mean(predict(fit, x) == y), 0.97)
})
