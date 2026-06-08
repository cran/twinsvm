test_that(".decision_grid returns expected columns and grid size", {
  set.seed(20)
  dat <- gen_moons(40, noise = 0.1)
  fit <- tsvm(dat$x, dat$y, kernel = "linear")
  grid <- .decision_grid(fit, dat$x, n = 17)

  expect_named(grid, c("x1", "x2", "class", "decision"))
  expect_equal(nrow(grid), 17 * 17)
  expect_s3_class(grid$class, "factor")
  expect_equal(levels(grid$class), levels(dat$y))
})

test_that("lift_plot returns ggplot objects for circles and moons", {
  set.seed(21)
  circles <- gen_circles(60, noise = 0.04)
  moons <- gen_moons(60, noise = 0.1)

  expect_s3_class(lift_plot(circles$x, circles$y, gamma = 1), "ggplot")
  expect_s3_class(lift_plot(moons$x, moons$y, gamma = 1), "ggplot")
})

test_that("lift_plotly errors cleanly when plotly is unavailable", {
  set.seed(22)
  dat <- gen_circles(30, noise = 0.04)
  old <- options(twinsvm.force_no_plotly = TRUE)
  on.exit(do.call(options, old), add = TRUE)

  expect_error(
    lift_plotly(dat$x, dat$y, gamma = 1),
    "requires the optional plotly package"
  )
})

test_that("compare_methods returns ggplots with all method labels", {
  set.seed(24)
  moons <- gen_moons(50, noise = 0.1)
  circles <- gen_circles(50, noise = 0.04)

  p_moons <- compare_methods(moons$x, moons$y, gamma = 1, c1 = 0.2, c2 = 0.2)
  p_circles <- compare_methods(circles$x, circles$y, gamma = 1, c1 = 0.2, c2 = 0.2)

  expect_s3_class(p_moons, "ggplot")
  expect_s3_class(p_circles, "ggplot")
  labels <- levels(p_moons$data$method)
  expect_equal(length(labels), 3)
  expect_true(any(grepl("^LS-TWSVM", labels)))
  expect_true(any(grepl("^QP-TWSVM", labels)))
  expect_true(any(grepl("^C-SVC SVM", labels)))
})

test_that("plot methods still return ggplots after grid refactor", {
  set.seed(23)
  dat <- gen_moons(50, noise = 0.1)
  fit_svm <- svms(dat$x, dat$y, kernel = "linear")
  fit_tsvm <- tsvm(dat$x, dat$y, kernel = "linear")

  expect_s3_class(plot(fit_svm), "ggplot")
  expect_s3_class(plot(fit_tsvm), "ggplot")
})
