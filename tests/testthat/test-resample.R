test_that("resample_accuracy returns a well-formed summary and per-split data", {
  set.seed(101)
  dat <- gen_moons(80, noise = 0.12)
  rs <- resample_accuracy(dat$x, dat$y, model = "both", reps = 5,
                          kernel = "rbf", gamma = 1, seed = 1)

  expect_s3_class(rs, "resample_accuracy")
  expect_equal(rs$reps, 5L)
  expect_setequal(as.character(rs$summary$model), c("tsvm", "svms"))
  expect_named(rs$summary, c("model", "mean", "sd", "se", "n"))
  expect_true(all(rs$summary$n == 5L))

  # se = sd / sqrt(n), and accuracies are valid proportions
  expect_equal(rs$summary$se, rs$summary$sd / sqrt(rs$summary$n))
  expect_equal(nrow(rs$accuracy), 10L) # 2 models x 5 reps
  expect_true(all(rs$accuracy$accuracy >= 0 & rs$accuracy$accuracy <= 1))
})

test_that("a fixed seed makes the evaluation reproducible", {
  set.seed(7)
  dat <- gen_moons(80, noise = 0.12)
  rs1 <- resample_accuracy(dat$x, dat$y, model = "tsvm", reps = 6,
                           kernel = "linear", seed = 42)
  rs2 <- resample_accuracy(dat$x, dat$y, model = "tsvm", reps = 6,
                           kernel = "linear", seed = 42)
  expect_equal(rs1$accuracy$accuracy, rs2$accuracy$accuracy)
})

test_that("arguments are forwarded only to the fitter that accepts them", {
  set.seed(3)
  dat <- gen_moons(80, noise = 0.12)
  # cost is svms-only, c1/c2 are tsvm-only; passing all three with model
  # "both" must not error -- each goes only where it is valid.
  expect_no_error(
    resample_accuracy(dat$x, dat$y, model = "both", reps = 4,
                      kernel = "rbf", gamma = 1,
                      cost = 1, c1 = 1, c2 = 1, seed = 5)
  )
})

test_that("multiclass responses are supported", {
  set.seed(9)
  x3 <- rbind(
    matrix(rnorm(40, -2, 0.3), ncol = 2),
    cbind(rnorm(20, 2, 0.3), rnorm(20, -2, 0.3)),
    matrix(rnorm(40, 2, 0.3), ncol = 2)
  )
  y3 <- factor(rep(c("a", "b", "c"), each = 20))
  rs <- resample_accuracy(x3, y3, model = "tsvm", reps = 4,
                          kernel = "linear", seed = 2)
  expect_equal(nrow(rs$accuracy), 4L)
  expect_true(rs$summary$mean > 0.8) # well-separated blobs
})

test_that("invalid settings are rejected", {
  set.seed(1)
  dat <- gen_moons(40, noise = 0.1)
  expect_error(resample_accuracy(dat$x, dat$y, reps = 1), "at least 2")
  expect_error(resample_accuracy(dat$x, dat$y, prop = 1), "between 0 and 1")
  expect_error(resample_accuracy(dat$x, dat$y, prop = 0), "between 0 and 1")

  x_one <- matrix(rnorm(6), ncol = 2)
  y_one <- factor(c("a", "a", "b"))  # class b has a single observation
  expect_error(resample_accuracy(x_one, y_one), "at least two observations")
})

test_that("plot method returns a ggplot", {
  set.seed(11)
  dat <- gen_moons(80, noise = 0.12)
  rs <- resample_accuracy(dat$x, dat$y, model = "both", reps = 5,
                          kernel = "rbf", gamma = 1, seed = 1)
  expect_s3_class(plot(rs), "ggplot")
})
