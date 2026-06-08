#' Confusion Matrix and Accuracy
#'
#' Predicts on `x` and returns a confusion matrix plus overall accuracy. This
#' works for both binary and one-vs-one multiclass `tsvm` and `svms` fits.
#'
#' @param object A fitted `tsvm`, `tsvm_multiclass`, `svms`, or
#'   `svms_multiclass` object.
#' @param x Numeric matrix or data frame of predictors.
#' @param y True class labels.
#'
#' @return A list with `table`, a `table(truth, predicted)` confusion matrix,
#'   and `accuracy`, the overall classification accuracy.
#' @family multiclass
#' @export
#'
#' @examples
#' set.seed(50)
#' x <- rbind(
#'   matrix(rnorm(20, -2, 0.2), ncol = 2),
#'   matrix(rnorm(20, 0, 0.2), ncol = 2),
#'   matrix(rnorm(20, 2, 0.2), ncol = 2)
#' )
#' y <- factor(rep(c("a", "b", "c"), each = 10))
#' fit <- tsvm(x, y, kernel = "linear")
#' confusion(fit, x, y)
confusion <- function(object, x, y) {
  if (is.null(object$levels)) {
    stop("`object` must be a fitted twinsvm model.", call. = FALSE)
  }
  pred <- stats::predict(object, x)
  if (length(y) != length(pred)) {
    stop("`y` must have the same length as the number of predictions.", call. = FALSE)
  }
  truth <- factor(y, levels = object$levels)
  if (anyNA(truth)) {
    stop("`y` must contain only classes seen during training.", call. = FALSE)
  }
  predicted <- factor(pred, levels = object$levels)
  tab <- table(truth = truth, predicted = predicted)
  list(
    table = tab,
    accuracy = mean(predicted == truth)
  )
}
