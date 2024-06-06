#' Title Variance estimation of regression coefficients
#'
#' @param object  an object inheriting from class RMST_Joint.
#' @export

vcov.jointModel_RMST <-
  function (object, ...) {
    out <- try(solve(object$Hessian), silent = TRUE)
    vmat <- if (!inherits(out, "try-error"))
      structure(out, dimnames = dimnames(object$Hessian))
    else
      structure(ginv(object$Hessian), dimnames = dimnames(object$Hessian))
    (vmat + t(vmat)) / 2
  }
