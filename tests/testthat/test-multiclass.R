three_blobs <- function(n = 90, sd = 0.18) {
  n_each <- n / 3L
  x <- rbind(
    matrix(rnorm(2L * n_each, mean = -2, sd = sd), ncol = 2),
    cbind(rnorm(n_each, mean = 2, sd = sd), rnorm(n_each, mean = -2, sd = sd)),
    matrix(rnorm(2L * n_each, mean = 2, sd = sd), ncol = 2)
  )
  y <- factor(rep(c("alpha", "beta", "gamma"), each = n_each))
  list(x = x, y = y)
}

test_that("binary tsvm path remains binary and numerically stable", {
  ref <- readRDS("ref-binary-tsvm.rds")
  dat <- list(x = ref$x, y = ref$y)
  fit <- tsvm(
    dat$x,
    dat$y,
    method = "ls",
    kernel = "rbf",
    gamma = 1.5,
    c1 = 0.4,
    c2 = 0.7
  )

  expect_s3_class(fit, "tsvm")
  expect_false(inherits(fit, "tsvm_multiclass"))
  expect_equal(fit$levels, ref$fit$levels)
  expect_equal(fit$method, ref$fit$method)
  expect_equal(fit$kernel, ref$fit$kernel)
  expect_equal(fit$u1, ref$fit$u1, tolerance = 1e-6)
  expect_equal(fit$b1, ref$fit$b1, tolerance = 1e-6)
  expect_equal(fit$u2, ref$fit$u2, tolerance = 1e-6)
  expect_equal(fit$b2, ref$fit$b2, tolerance = 1e-6)
  expect_identical(predict(fit, ref$x), ref$pred)
  expect_equal(predict(fit, ref$x, decision.values = TRUE), ref$decision, tolerance = 1e-6)
})

test_that("tsvm fits and predicts three-class OVO models", {
  set.seed(42)
  dat <- three_blobs()
  fit <- tsvm(dat$x, dat$y, kernel = "linear")
  pred <- predict(fit, dat$x)
  votes <- predict(fit, dat$x, type = "votes")

  expect_s3_class(fit, "tsvm_multiclass")
  expect_equal(fit$n_models, 3)
  expect_equal(nrow(fit$pairs), 3)
  expect_s3_class(pred, "factor")
  expect_equal(length(pred), nrow(dat$x))
  expect_equal(levels(pred), levels(dat$y))
  expect_gt(mean(pred == dat$y), 0.9)
  expect_equal(dim(votes), c(nrow(dat$x), 3L))
  expect_equal(rowSums(votes), rep(fit$n_models, nrow(dat$x)))
})

test_that("tsvm multiclass prediction handles one-row newdata", {
  set.seed(43)
  dat <- three_blobs()
  fit <- tsvm(dat$x, dat$y, kernel = "linear")
  pred <- predict(fit, dat$x[1, , drop = FALSE])

  expect_s3_class(pred, "factor")
  expect_equal(length(pred), 1)
  expect_equal(levels(pred), levels(dat$y))
})

test_that("tsvm multiclass decision values are rejected clearly", {
  set.seed(44)
  dat <- three_blobs()
  fit <- tsvm(dat$x, dat$y, kernel = "linear")

  expect_error(
    predict(fit, dat$x, decision.values = TRUE),
    "decision.values"
  )
})

test_that("binary svms path remains binary and predictively stable", {
  ref <- readRDS("ref-binary-svms.rds")
  set.seed(102)
  invisible(gen_moons(70, noise = 0.08))
  fit <- svms(
    ref$x,
    ref$y,
    kernel = "rbf",
    gamma = 1.5,
    cost = 0.8
  )

  expect_s3_class(fit, "svms")
  expect_false(inherits(fit, "svms_multiclass"))
  expect_equal(fit$levels, ref$fit$levels)
  expect_equal(fit$kernel, ref$fit$kernel)
  expect_equal(fit$cost, ref$fit$cost)
  expect_lte(abs(fit$n_support - ref$fit$n_support), 2)
  expect_identical(predict(fit, ref$x), ref$pred)
  expect_equal(predict(fit, ref$x, decision.values = TRUE), ref$decision, tolerance = 1e-2)
})

test_that("svms fits and predicts three-class OVO models", {
  set.seed(45)
  dat <- three_blobs()
  fit <- svms(dat$x, dat$y, kernel = "linear")
  pred <- predict(fit, dat$x)
  votes <- predict(fit, dat$x, type = "votes")

  expect_s3_class(fit, "svms_multiclass")
  expect_equal(fit$n_models, 3)
  expect_s3_class(pred, "factor")
  expect_equal(length(pred), nrow(dat$x))
  expect_equal(levels(pred), levels(dat$y))
  expect_gt(mean(pred == dat$y), 0.9)
  expect_equal(dim(votes), c(nrow(dat$x), 3L))
  expect_equal(rowSums(votes), rep(fit$n_models, nrow(dat$x)))
})

test_that("confusion returns table and accuracy", {
  set.seed(46)
  dat <- three_blobs()
  fit <- tsvm(dat$x, dat$y, kernel = "linear")
  cm <- confusion(fit, dat$x, dat$y)

  expect_named(cm, c("table", "accuracy"))
  expect_s3_class(cm$table, "table")
  expect_equal(dim(cm$table), c(3L, 3L))
  expect_gt(cm$accuracy, 0.9)
})

test_that("multiclass predictions are reproducible across repeated calls", {
  set.seed(47)
  dat <- three_blobs()
  fit <- tsvm(dat$x, dat$y, kernel = "linear")

  expect_identical(predict(fit, dat$x), predict(fit, dat$x))
  expect_identical(
    predict(fit, dat$x, type = "votes"),
    predict(fit, dat$x, type = "votes")
  )
})

test_that("degenerate OVO pairs fail with informative pair message", {
  x <- rbind(
    matrix(rnorm(8, mean = -2, sd = 0.1), ncol = 2),
    matrix(rnorm(8, mean = 2, sd = 0.1), ncol = 2),
    c(0, 0)
  )
  y <- factor(c(rep("a", 4), rep("b", 4), "c"))

  expect_error(
    tsvm(x, y, kernel = "linear"),
    "pair `a` vs `c`"
  )
})

test_that("empty levels are dropped before OVO pair construction", {
  set.seed(48)
  dat <- three_blobs()
  y <- factor(dat$y, levels = c(levels(dat$y), "empty"))
  fit <- tsvm(dat$x, y, kernel = "linear")

  expect_s3_class(fit, "tsvm_multiclass")
  expect_equal(fit$levels, levels(dat$y))
  expect_equal(fit$n_models, 3)
})

test_that("one-class responses and column mismatch error clearly", {
  set.seed(49)
  dat <- three_blobs()
  expect_error(
    tsvm(dat$x, factor(rep("only", nrow(dat$x))), kernel = "linear"),
    "at least two classes"
  )

  fit <- tsvm(dat$x, dat$y, kernel = "linear")
  expect_error(
    predict(fit, cbind(dat$x, 1)),
    "must have 2 columns"
  )
})
