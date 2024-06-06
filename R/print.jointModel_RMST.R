#' Title output result
#'
#' @param out an object inheriting from class RMST_Joint.
#' @export

print.jointModel_RMST<-function(out,...){
  print(out$coefficients,...)
}
