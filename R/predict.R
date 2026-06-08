#' Predict from a Standard SVM
#'
#' @param object A fitted `svms` object.
#' @param newdata Numeric matrix or data frame.
#' @param decision.values If `TRUE`, return raw decision values instead of
#'   class labels.
#' @param ... Unused.
#'
#' @return A factor of predicted classes, or a numeric vector when
#'   `decision.values = TRUE`.
#' @export
#'
#' @examples
#' set.seed(2)
#' dat <- gen_moons(30)
#' fit <- svms(dat$x, dat$y)
#' predict(fit, dat$x, decision.values = TRUE)
predict.svms <- function(object, newdata, decision.values = FALSE, ...) {
  x <- check_newdata(object, newdata)
  if (length(object$support_alpha) == 0L) {
    decision <- rep(object$b, nrow(x))
  } else {
    k <- kernel_matrix_cpp(
      x, object$support_vectors, object$kernel, object$gamma,
      object$degree, object$coef0
    )
    decision <- drop(k %*% (object$support_alpha * object$support_y) + object$b)
  }
  if (decision.values) {
    return(decision)
  }
  label_from_decision(decision, object$levels)
}
