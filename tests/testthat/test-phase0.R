test_that("data generators return valid two-class data", {
  set.seed(1)
  dat <- gen_moons(25, noise = 0.1)
  expect_equal(nrow(dat$x), 25)
  expect_equal(ncol(dat$x), 2)
  expect_s3_class(dat$y, "factor")
  expect_equal(levels(dat$y), c("B", "A"))
})

test_that("svms fits and predicts with stable shape", {
  set.seed(2)
  dat <- gen_moons(60, noise = 0.05)
  fit <- svms(dat$x, dat$y, kernel = "linear", cost = 1)
  pred <- predict(fit, dat$x)
  dec <- predict(fit, dat$x, decision.values = TRUE)
  expect_s3_class(fit, "svms")
  expect_equal(length(pred), nrow(dat$x))
  expect_equal(levels(pred), levels(dat$y))
  expect_type(dec, "double")
  expect_equal(length(dec), nrow(dat$x))
})

test_that("linear LSTSVM class planes prefer their own class on blobs", {
  set.seed(3)
  x_b <- matrix(rnorm(80, mean = -2, sd = 0.25), ncol = 2)
  x_a <- matrix(rnorm(80, mean = 2, sd = 0.25), ncol = 2)
  x <- rbind(x_b, x_a)
  y <- factor(c(rep("B", 40), rep("A", 40)), levels = c("B", "A"))
  fit <- tsvm(x, y)
  pred <- predict(fit, x)
  expect_s3_class(fit, "tsvm")
  expect_gt(mean(pred == y), 0.97)

  f1 <- abs(drop(x %*% fit$w1 + fit$b1)) / fit$norm1
  f2 <- abs(drop(x %*% fit$w2 + fit$b2)) / fit$norm2
  expect_lt(mean(f1[y == "A"]), mean(f1[y == "B"]))
  expect_lt(mean(f2[y == "B"]), mean(f2[y == "A"]))
})

test_that("kernel LSTSVM improves over linear LSTSVM on moons", {
  set.seed(4)
  dat <- gen_moons(120, noise = 0.12)
  linear <- tsvm(dat$x, dat$y, kernel = "linear")
  kernel <- tsvm(dat$x, dat$y, kernel = "rbf", gamma = 2, c1 = 0.1, c2 = 0.1)
  acc_linear <- mean(predict(linear, dat$x) == dat$y)
  acc_kernel <- mean(predict(kernel, dat$x) == dat$y)
  expect_gt(acc_kernel, acc_linear + 0.05)
  expect_gt(acc_kernel, 0.9)
  expect_equal(length(predict(kernel, dat$x, decision.values = TRUE)), nrow(dat$x))
})

test_that("original QP twin SVM is close to least-squares twin SVM", {
  set.seed(5)
  dat <- gen_moons(90, noise = 0.1)
  fit_ls <- tsvm(dat$x, dat$y, method = "ls", kernel = "rbf", gamma = 2, c1 = 0.2, c2 = 0.2)
  fit_qp <- tsvm(dat$x, dat$y, method = "twin", kernel = "rbf", gamma = 2, c1 = 0.2, c2 = 0.2)
  acc_ls <- mean(predict(fit_ls, dat$x) == dat$y)
  acc_qp <- mean(predict(fit_qp, dat$x) == dat$y)
  expect_gt(acc_qp, 0.85)
  expect_lt(abs(acc_ls - acc_qp), 0.2)

  x_b <- matrix(rnorm(60, mean = -1.5, sd = 0.2), ncol = 2)
  x_a <- matrix(rnorm(60, mean = 1.5, sd = 0.2), ncol = 2)
  x <- rbind(x_b, x_a)
  y <- factor(c(rep("B", 30), rep("A", 30)), levels = c("B", "A"))
  fit_linear <- tsvm(x, y, method = "twin", kernel = "linear")
  expect_gt(mean(predict(fit_linear, x) == y), 0.95)
})

test_that("cv_tsvm returns complete grid and best parameters", {
  set.seed(6)
  dat <- gen_moons(60, noise = 0.1)
  cv <- cv_tsvm(
    dat$x, dat$y,
    c1_grid = c(0.1, 1),
    c2_grid = c(0.1, 1),
    gamma_grid = c(1, 2),
    k = 3,
    kernel = "rbf"
  )
  expect_s3_class(cv, "cv_tsvm")
  expect_equal(nrow(cv$results), 8)
  expect_true(cv$best_params$c1 %in% c(0.1, 1))
  expect_true(cv$best_params$c2 %in% c(0.1, 1))
  expect_true(cv$best_params$gamma %in% c(1, 2))
  expect_s3_class(plot(cv), "ggplot")
})

test_that("visualization functions return plot objects", {
  set.seed(7)
  dat <- gen_moons(50, noise = 0.1)
  fit_svm <- svms(dat$x, dat$y, kernel = "linear")
  fit_tsvm <- tsvm(dat$x, dat$y, kernel = "linear")
  expect_s3_class(plot(fit_svm), "ggplot")
  expect_s3_class(plot(fit_tsvm), "ggplot")
  expect_s3_class(kernel_lift(dat$x, dat$y), "ggplot")
  old_wd <- setwd(tempdir())
  on.exit({
    setwd(old_wd)
    unlink(file.path(tempdir(), "gganim_plot*.png"))
  }, add = TRUE)
  anim <- morph_boundary(
    dat$x, dat$y,
    param = "gamma",
    range = c(0.5, 1),
    model = "tsvm",
    n = 3,
    kernel = "rbf"
  )
  expect_true(inherits(anim, "gganim"))
})
