#' Accuracy and Its Spread over Repeated Random Splits
#'
#' Evaluates a classifier by repeated random holdout: the data are split into
#' training and test sets `reps` times, the model is refitted on each training
#' set, and test accuracy is recorded. Reporting the mean **together with its
#' between-split standard deviation and standard error** turns a single accuracy
#' number into an estimate with a stated uncertainty, so that differences
#' between methods can be read against the noise rather than taken at face
#' value.
#'
#' The splits are drawn once, up front, and shared by every model, so when
#' `model = "both"` the twin SVM and the standard SVM are compared on identical
#' training and test sets (a paired comparison). Splits are stratified by class,
#' so each class keeps roughly its overall proportion in both parts.
#'
#' @param x Numeric matrix or data frame of predictors.
#' @param y Response factor with at least two classes.
#' @param model Which classifier(s) to evaluate: `"both"` (the default,
#'   [tsvm()] and [svms()] on the same splits), `"tsvm"`, or `"svms"`.
#' @param reps Number of random splits (at least 2; the spread is estimated
#'   across them).
#' @param prop Proportion of observations used for training in each split,
#'   between 0 and 1.
#' @param seed Optional integer seed. Set it for reproducible splits; because
#'   the splits are fixed before any model is fitted, the comparison is stable
#'   even though the SMO solver in [svms()] itself draws random numbers.
#' @param ... Further arguments passed to the fitter(s), e.g. `kernel`,
#'   `gamma`, `cost` (standard SVM), or `c1`, `c2`, `method` (twin SVM). Each
#'   argument is forwarded only to the fitter(s) that accept it, so a single
#'   call such as `kernel = "rbf", gamma = 1` works for `model = "both"`.
#'
#' @return An object of class `resample_accuracy`, a list with
#'   \describe{
#'     \item{`summary`}{a data frame with one row per model giving `mean`, `sd`,
#'       `se` and `n` of the test accuracy;}
#'     \item{`accuracy`}{a long data frame of the per-split accuracies
#'       (`model`, `rep`, `accuracy`);}
#'     \item{`reps`, `prop`, `models`, `call`}{the settings used.}
#'   }
#'   The [plot()][plot.resample_accuracy] method draws the accuracy
#'   distributions as boxplots.
#'
#' @importFrom stats sd predict
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- gen_moons(120, noise = 0.15)
#' rs <- resample_accuracy(dat$x, dat$y, model = "both",
#'                         reps = 10, kernel = "rbf", gamma = 1, seed = 1)
#' rs                 # mean, sd and se per model
#' plot(rs)           # accuracy spread as boxplots
resample_accuracy <- function(x, y, model = c("both", "tsvm", "svms"),
                              reps = 50L, prop = 0.7, seed = NULL, ...) {
  x <- as_numeric_matrix(x)
  y <- check_class_factor(y, nrow(x))
  model <- match.arg(model)
  models <- if (model == "both") c("tsvm", "svms") else model

  if (!is.numeric(reps) || length(reps) != 1L || is.na(reps) || reps < 2L) {
    stop("`reps` must be a single integer of at least 2.", call. = FALSE)
  }
  reps <- as.integer(reps)
  if (!is.numeric(prop) || length(prop) != 1L || is.na(prop) ||
      prop <= 0 || prop >= 1) {
    stop("`prop` must be a single number strictly between 0 and 1.",
         call. = FALSE)
  }
  if (min(table(y)) < 2L) {
    stop("each class must have at least two observations to be split.",
         call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Draw all splits up front so every model sees identical train/test sets.
  splits <- replicate(reps, stratified_train_idx(y, prop), simplify = FALSE)

  extra <- list(...)
  n <- nrow(x)
  all_idx <- seq_len(n)

  per_model <- lapply(models, function(m) {
    fitter <- if (m == "tsvm") tsvm else svms
    args_ok <- extra[names(extra) %in% names(formals(fitter))]
    acc <- vapply(splits, function(tr_idx) {
      te_idx <- all_idx[-tr_idx]
      fit <- do.call(fitter, c(
        list(x = x[tr_idx, , drop = FALSE], y = y[tr_idx]), args_ok
      ))
      mean(predict(fit, x[te_idx, , drop = FALSE]) == y[te_idx])
    }, numeric(1L))
    data.frame(model = m, rep = seq_len(reps), accuracy = acc,
               stringsAsFactors = FALSE)
  })

  acc_df <- do.call(rbind, per_model)
  acc_df$model <- factor(acc_df$model, levels = models)

  summary_df <- do.call(rbind, lapply(models, function(m) {
    a <- acc_df$accuracy[acc_df$model == m]
    s <- sd(a)
    data.frame(model = m, mean = mean(a), sd = s,
               se = s / sqrt(length(a)), n = length(a),
               stringsAsFactors = FALSE)
  }))
  rownames(summary_df) <- NULL

  structure(
    list(
      summary = summary_df,
      accuracy = acc_df,
      reps = reps,
      prop = prop,
      models = models,
      call = match.call()
    ),
    class = "resample_accuracy"
  )
}

# Stratified training indices: sample a `prop` fraction within each class,
# keeping at least one observation of that class in each of train and test.
stratified_train_idx <- function(y, prop) {
  by_class <- split(seq_along(y), y)
  picked <- lapply(by_class, function(ii) {
    k <- floor(prop * length(ii))
    k <- min(max(k, 1L), length(ii) - 1L)
    sample(ii, k)
  })
  sort(unlist(picked, use.names = FALSE))
}

#' @export
print.resample_accuracy <- function(x, ...) {
  cat("Repeated random holdout (", x$reps, " splits, ",
      format(x$prop), " train)\n", sep = "")
  s <- x$summary
  out <- data.frame(
    model = s$model,
    mean = formatC(s$mean, digits = 4, format = "f"),
    sd = formatC(s$sd, digits = 4, format = "f"),
    se = formatC(s$se, digits = 4, format = "f"),
    stringsAsFactors = FALSE
  )
  print(out, row.names = FALSE)
  invisible(x)
}

#' Plot the Accuracy Spread from Repeated Splits
#'
#' Draws the per-split test accuracies as a boxplot, one box per model, so the
#' dispersion is visible alongside the mean.
#'
#' @param x A `resample_accuracy` object.
#' @param ... Unused.
#'
#' @return A `ggplot` object.
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' set.seed(2)
#' dat <- gen_moons(120, noise = 0.15)
#' rs <- resample_accuracy(dat$x, dat$y, reps = 10, kernel = "rbf", gamma = 1)
#' plot(rs)
plot.resample_accuracy <- function(x, ...) {
  ggplot2::ggplot(
    x$accuracy,
    ggplot2::aes(.data$model, .data$accuracy, fill = .data$model)
  ) +
    ggplot2::geom_boxplot(width = 0.55, alpha = 0.75,
                          outlier.size = 0.8) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 23,
                          size = 2.4, fill = "white") +
    ggplot2::labs(x = NULL, y = "Test accuracy", fill = NULL) +
    ggplot2::guides(fill = "none") +
    .twin_theme()
}
