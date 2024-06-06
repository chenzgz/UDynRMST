#' Title Summary of output results
#'
#' @param object an object inheriting from class RMST_Joint.
#' @export

summary.jointModel_RMST<-function(object){
  VarCov <- vcov(object)
  betas <- object$coefficients$betas
  indY <- seq(1, length(betas))
  sds <- sqrt(diag(VarCov[indY, indY]))
  coefsY <- cbind("Value" = betas, "Std.Err" = sds, "z-value" = betas / sds, "Lower_CI"=betas-1.96*sds,"Upper_CI"=betas+1.96*sds,
                  "p-value" = 2 * pnorm(abs(betas / sds), lower.tail = FALSE))
  gammas <- c(object$coefficients$gammas,
              "Assoct" = as.vector(object$coefficients$alpha),
              "Assoct.s" = as.vector(object$coefficients$Dalpha))
  indT <- grep("T.", colnames(VarCov), fixed = TRUE)
  indT<-head(indT,-1)
  jj <- grep("Assoct[!^\\.s]", names(gammas))
  ii <- setdiff(grep("Assoct", names(gammas)), jj)
  if (length(ii) > 1) {
    nn <- names(object$coefficients$alpha)
    names(gammas)[ii] <- if (length(nn) == 1) "Assoct" else {
      if (nn[1] == "")
        c("Assoct", paste("Assoct", nn[-1], sep = ":"))
      else
        paste("Assoct", nn, sep = ":")
    }
  }
  if (length(jj) > 1) {
    nn <- names(object$coefficients$Dalpha)
    names(gammas)[jj] <- if (length(nn) == 1) "Assoct.s" else {
      if (nn[1] == "")
        c("Assoct.s", paste("Assoct.s", nn[-1], sep = ":"))
      else
        paste("Assoct.s", nn, sep = ":")
    }
  }
  sds <- if (length(indT) > 1) sqrt(diag(VarCov[indT, indT])) else sqrt(VarCov[indT, indT])
  if (!is.null(object$scaleWB))
    sds <- c(sds, NA)
  coefsT <- cbind("Value" = gammas, "Std.Err" = sds, "z-value" = gammas / sds,"Lower_CI"=gammas-1.96*sds,"Upper_CI"=gammas+1.96*sds,
                  "p-value" = 2 * pnorm(abs(gammas / sds), lower.tail = FALSE))
  out1<-list("CoefTable-Long" = coefsY, "CoefTable-Event" = coefsT,"sigma"=object$coefficients$sigma,"delta"=object$coefficients$delta,"D"=object$coefficients$D)
  class(out1) <- "summary.jointModel_RMST"
  out1
}
