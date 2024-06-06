#' Title Dynamic prediction in joint model based on restricted mean survival time (RMST)
#' @description
#' This function computes the the survival time after the last observed time for which a longitudinal measurement was available.
#' @usage survJM(object,newdata,idVar, simulate = TRUE, survTimes = NULL, last.time = NULL, M = 200, CI.levels = c(0.025, 0.975), scale = 1.6)
#' @param object an object inheriting from class RMST_Joint.
#' @param newdata a data frame that contains the longitudinal and covariate information for the subjects for which prediction of dynamic RMST is required.
#' @param idVar the name of the variable in newdata that identifies the different subjects.
#' @param simulate logical; if TRUE, a Monte Carlo approach is used to estimate dynamic RMST. If FALSE, a first order estimator is used instead.
#' @param survTimes a numeric vector of times for which prediction dynamic RMST are to be computed
#' @param last.time a numeric vector. This specifies the known time at which each of the subjects in newdat was known to be alive. If NULL, then this is automatically taken as the last time each subject provided a longitudinal measurement.
#' @param M integer denoting how many Monte Carlo samples to use
#' @param CI.levels a numeric vector of length two that specifies which quantiles to use for the calculation of confidence interval for the predicted dynamic RMST.
#' @param scale a numeric scalar that controls the acceptance rate of the Metropolis-Hastings algorithm
#'
#' @return A list of class survJM with components:
#' \describe{
#'  \item{summaries}{a list with elements numeric matrices with numeric summaries of the predicted survival time for each subject.}
#'  \item{survTimes}{a copy of the survTimes argument.}
#'  \item{last.time}{a numeric vector with the time of the last available longitudinal measurement of each subject.}
#'  \item{obs.times}{a list with elements numeric vectors denoting the timings of the longitudinal measurements for each subject.}
#'  \item{y}{a list with elements numeric vectors denoting the longitudinal responses for each subject.}
#'  \item{fitted.times}{a list with elements numeric vectors denoting the timings of the longitudinal measurements for each subject and last.time.}
#'  \item{fitted.y}{Measures of longitudinal covariates from longitudinal model fits}
#'  \item {ry} {Range of values for the measurement of longitudinal covariates}}
#' @export
#' @examples
#'library(UDynRMST)
#'data(pbc2)
#'data(pbc2.id)
#'head(pbc2)
#'head(pbc2.id)
#'pbc2$status<-as.numeric(pbc2$status=="dead")
#'pbc2$drug<-as.numeric(pbc2$drug=="D-penicil")
#'pbc2$sex<-as.numeric(pbc2$sex=="female")
#'pbc2.id$status<-as.numeric(pbc2.id$status=="dead")
#'pbc2.id$drug<-as.numeric(pbc2.id$drug=="D-penicil")
#'pbc2.id$sex<-as.numeric(pbc2.id$sex=="female")
#'tau=max(pbc2[which(pbc2$status==1),]$years)
#'pbc2$time<-pmin(pbc2$years,tau)
#'pbc2$status1<-ifelse(pbc2$time==pbc2$years,pbc2$status,0)
#'indd<-pbc2$year>tau
#'data<-pbc2[!indd,]
#'pbc2.id$time<-pmin(pbc2.id$years,tau)
#'pbc2.id$status1<-ifelse(pbc2.id$time==pbc2.id$years,pbc2.id$status,0)
#'data.id<-pbc2.id
#'data<-data[order(data$id),]
#'data.id<-data.id[order(data.id$id),]
#'fit_nonlinear3 <- lme(log(serBilir) ~ ns(year,3)+ age + sex + drug,
#'random =list(id = pdDiag(form = ~ ns(year, 3))), data =data)
#'p<-pseudomean(data.id$time,data.id$status1,tau)
#'fitSURV<-geeglm(p~drug+sex+age,data=data.id,id= id)
#'object5<-jointModelRMST(fit_nonlinear3,fitSURV,data.id$time,"year")
#'newdata_5 <- pbc2[pbc2$id ==5, ]
#'#s<-survJM(object = object5, newdata = newdata_5,idVar = "id",simulate =FALSE)
#'survRes_5 <- list()
#'for (o in 1:dim(newdata_5)[1]) {
#'  newdataDyn_5 <- newdata_5[1:o, ]
#'  survRes_5[[o]] <-  survJM(object = object5, newdata = newdataDyn_5,idVar = "id",
#'                            simulate = TRUE, survTimes = NULL, last.time = NULL,
#'                            M =1000, CI.levels = c(0.025, 0.975), scale = 1.6)
#'}

#'col3="darkorange2"
#'col1="steelblue4"
#'col2="lightcoral"
#'i=6
#'par(par(mar = c(5,4,4,2.1)),oma=c(0.5,1,0.5,1))
#'rng <- range(survRes_5[[i]]$obs.times[[1]], survRes_5[[i]]$survTimes)
#'plot(survRes_5[[i]]$obs.times[[1]],survRes_5[[i]]$y[[1]],xlim=rng,ylim=survRes_5[[i]]$ry,pch=8,ylab="",xlab="",
#'     xaxt="n",yaxs = "i",yaxt="n",lwd=3,col=col1)
#'axis(1,0:15,cex.axis=1.5,lwd=1.6,font.axis=2,tcl=-0.4,hadj=0.9)
#'axis(2,cex.axis=1.5,lwd=1.6,font.axis=2,tcl=-0.4,hadj=0.9,col.axis=col1,col.ticks =col1)
#'title(ylab="log serum bilirubin",line=2.2, cex.lab=1.5,las=0,font.lab=2,col.lab=col1)
#'title(xlab="Prediction time in years (s)",line=2.2, cex.lab=1.5,las=0,font.lab=2)
#'lines(survRes_5[[i]]$fitted.times[[1]],survRes_5[[i]]$fitted.y[[1]],lwd=3)
#'abline(v = survRes_5[[i]]$last.time[1], lty = 3,lwd=3)
#'par(new = TRUE)
#'r.=survRes_5[[i]]$summaries[[1]]
#'matplot(r.[, 1], r.[, -c(1,3), drop = FALSE], type = "l", col =  c(col3,col3, col3), lwd = c(6, 6, 6),
#'        lty = c(1, 3, 3),xlim = rng,ylim=c(0,20),
#'        ylab = "", xlab = "", axes = FALSE, yaxs = "i")
#'axis(4, las = 2, cex.axis = 1.5,line=0,lwd=1.6,font.axis=2,tcl=-0.3,hadj=0.6,col.axis=col3,col.ticks =col3)
#'mtext("Dynamic RMST",side=4,line=1.8, cex=1.5,las=0,font=2,col=col3)
#'title(list("c                                   ",font=7,cex=2.5))
#'title(list(paste0("Patient A ", ifelse(i==1,paste(" ","1st"),ifelse(i==2,paste(" ","2nd"),paste("  ",i,"th")))),font=7 ,cex=2))

survJM<-function(object,newdata,idVar, simulate = TRUE, survTimes = NULL, last.time = NULL, M = 200, CI.levels = c(0.025, 0.975), scale = 1.6){
  if (is.null(survTimes) || !is.numeric(survTimes))
    survTimes <- seq(min((object$yy$Time)),
                     max((object$yy$Time)) + 0.1, length.out = 35)
  timeVar <- object$timeVar
  TermsX <- object$termsYx
  TermsZ <- object$termsYz
  LongFormat <- FALSE
  mfX <- model.frame(TermsX, data = newdata)
  mfZ <- model.frame(TermsZ, data = newdata)
  formYx <- reformulate(attr(delete.response(TermsX), "term.labels"))
  formYz <- object$formYz
  na.ind <- as.vector(attr(mfX, "na.action"))
  na.ind <- if (is.null(na.ind)) {
    rep(TRUE, nrow(newdata))
  } else {
    !seq_len(nrow(newdata)) %in% na.ind
  }
  id <- as.numeric(unclass(newdata[[idVar]]))
  id <- id. <- match(id, unique(id))
  id <- id[na.ind]
  y <- model.response(mfX)
  X <- model.matrix(formYx, mfX)
  Z <- model.matrix(formYz, mfZ)[na.ind, , drop = FALSE]
  TermsT <- object$termsT
  data.id <-newdata[tapply(row.names(newdata), id, tail, n = 1L),]
  idT <- data.id[[idVar]]
  idT <- match(idT, unique(idT))
  mfT <- model.frame(delete.response(TermsT), data = data.id)
  formT <- if (!is.null(kk <- attr(TermsT, "specials")$strata)) {
    strt <- eval(attr(TermsT, "variables"), data.id)[[kk]]
    tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
    reformulate(attr(tt, "term.labels"))
  } else if (!is.null(kk <- attr(TermsT, "specials")$cluster)) {
    tt <- drop.terms(TermsT, kk - 1, keep.response = FALSE)
    reformulate(attr(tt, "term.labels"))
  } else {
    tt <- attr(delete.response(TermsT), "term.labels")
    if (length(tt)) reformulate(tt) else reformulate("1")
  }
  W <- model.matrix(formT, mfT)
  obs.times <- split(newdata[[timeVar]][na.ind], id)
  last.time <- if (is.null(last.time)) {
    tapply(newdata[[timeVar]], id., tail, n = 1)
  } else if (is.character(last.time) && length(last.time) == 1) {
    tapply(newdata[[last.time]], id., tail, n = 1)
  } else if (is.numeric(last.time) && length(last.time) == nrow(data.id)) {
    last.time
  } else {
    stop("\nnot appropriate value for 'last.time' argument.")
  }

  times.to.pred <- lapply(last.time, function (t) survTimes[survTimes > t])
  n <- object$n
  n.tp <- length(last.time)
  ncx <- ncol(X)
  ncz <- ncol(Z)
  ncww <- ncol(W)
  betas <- object$coefficients[['betas']]
  sigma <- object$coefficients[['sigma']]
  D <- object$coefficients[['D']]
  diag.D <- ncol(D) == 1 & nrow(D) > 1
  D <- if (diag.D) diag(c(D)) else D
  gammas <- object$coefficients[['gammas']]
  alpha <- object$coefficients[['alpha']]
  delta<-object$coefficients[['delta']]
  list.thetas <- list(betas = betas, log.sigma = log(sigma), gammas = gammas, alpha = alpha,
                      delta =log( delta),   D = if (diag.D) log(diag(D)) else JM:::chol.transf(D))
  list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
  thetas <- unlist(as.relistable(list.thetas))
  Var.thetas <- vcov(object)
  environment(log.posterior.b) <- environment(S.b) <- environment()
  # construct model matrices to calculate the survival functions
  #obs.times.surv <- split(data.id[[timeVar]], idT)
  #survMats <- survMats.last <- vector("list", n.tp)

  # calculate the Empirical Bayes estimates and their (scaled) variance
  modes.b <- matrix(0, n.tp, ncz)
  Vars.b <- vector("list", n.tp)
  for (i in seq_len(n.tp)) {
    betas.new <- betas
    sigma.new <- sigma
    D.new <- D
    gammas.new <- gammas
    alpha.new <- alpha
    delta.new <- delta
    ff <- function (b, y, tt, i) -log.posterior.b(b, y, time=tt, ii = i)
    opt <- try(optim(rep(0, ncz), ff, y = y, tt=last.time[i], i = i,
                     method = "BFGS", hessian = TRUE), TRUE)
    if (inherits(opt, "try-error")) {
      gg <- function (b, y, tt, i) JM:::cd(b, ff, y = y, tt=tt, i = i)
      opt <- optim(rep(0, ncz), ff, gg, y = y, tt=last.time[i],
                   i = i, method = "BFGS", hessian = TRUE)
    }
    modes.b[i, ] <- opt$par
    Vars.b[[i]] <- scale * solve(opt$hessian)
  }
  if (!simulate) {
    res <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
      S.last <- S.b(last.time[i], modes.b[i, ], i)
      S.pred <- numeric(length(times.to.pred[[i]]))
      for (l in seq_along(S.pred))
        S.pred[l] <- S.b(times.to.pred[[i]][l], modes.b[i, ], i)
      res[[i]] <- cbind(times = c(last.time[i],times.to.pred[[i]]), predSurv = c(S.last-last.time[i],S.pred-times.to.pred[[i]]))
      rownames(res[[i]]) <- seq_along(c(S.last,S.pred))
    }
  } else {
    out <- vector("list", M)
    success.rate <- matrix(FALSE, M, n.tp)
    b.old <- b.new <- modes.b
    if (n.tp == 1)
      dim(b.old) <- dim(b.new) <- c(1, ncz)
    for (m in 1:M) {
      # Step 1: simulate new parameter value
      set.seed(0104+m)
      thetas.new <- mvrnorm(1, thetas, Var.thetas)
      thetas.new <- relist(thetas.new, skeleton = list.thetas)
      betas.new <- thetas.new$betas
      sigma.new <- exp(thetas.new$log.sigma)
      gammas.new <- thetas.new$gammas
      alpha.new <- thetas.new$alpha
      delta.new<-exp(thetas.new$delta)
      D.new <- thetas.new$D
      D.new <- if (diag.D) exp(D.new) else JM:::chol.transf(D.new)
      SS <- vector("list", n.tp)
      for (i in seq_len(n.tp)) {
        # Step 2: simulate new random effects values
        proposed.b <- JM:::rmvt(1, modes.b[i, ], Vars.b[[i]], 4)
        dmvt.old <- JM:::dmvt(b.old[i, ], modes.b[i, ], Vars.b[[i]], 4, TRUE)
        dmvt.proposed <- JM:::dmvt(proposed.b, modes.b[i, ], Vars.b[[i]], 4, TRUE)
        a <- min(exp(log.posterior.b(proposed.b, y, last.time[[i]], ii = i) + dmvt.old -
                       log.posterior.b(b.old[i, ], y, last.time[[i]], ii = i) - dmvt.proposed), 1)
        ind <- runif(1) <= a
        success.rate[m, i] <- ind
        if (!is.na(ind) && ind)
          b.new[i, ] <- proposed.b
        # Step 3: compute Pr(T > t_k | T > t_{k - 1}; theta.new, b.new)

        S.last <- S.b(last.time[i], b.new[i, ], i)
        S.pred <- numeric(length(times.to.pred[[i]]))
        for (l in seq_along(S.pred))
          S.pred[l] <- S.b(times.to.pred[[i]][l], b.new[i, ], i)
        SS[[i]]<-c(S.last-last.time[i],S.pred-times.to.pred[[i]])


      }
      b.old <- b.new
      out[[m]] <- SS
    }
    res <- vector("list", n.tp)
    for (i in seq_len(n.tp)) {
      rr <- sapply(out, "[[", i)
      if (!is.matrix(rr))
        rr <- rbind(rr)
      res[[i]] <- cbind(
        times =c(last.time[i],times.to.pred[[i]]),
        "Mean" = rowMeans(rr, na.rm = TRUE),
        "Median" = apply(rr, 1, median, na.rm = TRUE),
        "Lower" = apply(rr, 1, quantile, probs = CI.levels[1], na.rm = TRUE),
        "Upper" = apply(rr, 1, quantile, probs = CI.levels[2], na.rm = TRUE)
      )
      rownames(res[[i]]) <- as.character(seq_len(NROW(res[[i]])))
    }
  }
  y <- split(y, id)
  newdata. <- do.call(rbind, mapply(function (d, t) {
    d. <- rbind(d, d[nrow(d), ])
    d.[[timeVar]][nrow(d.)] <- t
    d.
  }, split(newdata, id.), last.time, SIMPLIFY = FALSE))
  id. <- as.numeric(unclass(newdata.[[idVar]]))
  id. <- match(id., unique(id.))
  mfX. <- model.frame(delete.response(TermsX), data = newdata.)
  mfZ. <- model.frame(TermsZ, data = newdata.)
  X. <- model.matrix(formYx, mfX.)
  Z. <- model.matrix(formYz, mfZ.)
  fitted.y <- split(c(X. %*% betas) + rowSums(Z. * modes.b[id., , drop = FALSE]), id.)
  names(res) <-names(y) <- names(last.time) <- names(obs.times) <- unique(unclass(newdata[[idVar]]))
  res <- list(summaries = res, survTimes = survTimes, last.time = last.time,
              obs.times = obs.times, y =y,
              fitted.times = split(newdata.[[timeVar]], factor(newdata.[[idVar]])),
              fitted.y = fitted.y, ry = range(object$y$y, na.rm = TRUE))
  if (simulate) {
    res$full.results <- out
    res$success.rate <- success.rate
    rm(list = ".Random.seed", envir = globalenv())
  }
  class(res) <- "survJM"
  res
}
S.b <-
  function (t, b, ii) {
    idT.i <- idT %in% ii
    eta.tw <- if (!is.null(W))
      as.vector(W[ii, , drop = FALSE] %*% gammas.new)
    else 0
    data.id[[timeVar]]<- t
    mfX.id<-model.frame(TermsX,data.id)
    mfZ.id<-model.frame(TermsZ,data.id)
    Xtime.i <- model.matrix(formYx, mfX.id)[idT.i, , drop = FALSE]
    Ztime.i <- model.matrix(formYz, mfZ.id)[idT.i, , drop = FALSE]
    Y<-as.vector(Xtime.i %*% betas.new+rowSums(Ztime.i*rep(b,each=nrow(Ztime.i))))
    mu.T <- eta.tw+alpha.new*Y
    mu.T
  }

#b=rep(0, ncz)
#y=y
#time=last.time[84]
#ii=84

log.posterior.b <-
  function (b, y, time, ii) {
    id.i <- id %in% ii
    idT.i <- idT %in% ii
    X.i <- X[id.i, , drop = FALSE]
    Z.i <- Z[id.i, , drop = FALSE]
    data.id[[timeVar]]<- time
    mfX.id<-model.frame(TermsX,data.id)
    mfZ.id<-model.frame(TermsZ,data.id)
    Xtime.i <- model.matrix(formYx, mfX.id)[idT.i, , drop = FALSE]
    Ztime.i <- model.matrix(formYz, mfZ.id)[idT.i, , drop = FALSE]
    mu.y <- as.vector(X.i %*% betas.new) + rowSums(Z.i * rep(b, each = nrow(Z.i)))
    logNorm <- dnorm(y[id.i], mu.y, sigma.new, TRUE)
    log.p.yb <- sum(logNorm)
    log.p.b <- JM:::dmvnorm(b, rep(0, ncol(Z)), D.new, TRUE)
    eta.tw <- if (!is.null(W)) {
      if (!LongFormat)
        as.vector(W[ii, , drop = FALSE] %*% gammas.new)
      else
        as.vector(W[idT.i, , drop = FALSE] %*% gammas.new)
    } else 0
    Y<-as.vector( Xtime.i %*% betas.new+rowSums(Ztime.i*rep(b,each=nrow(Ztime.i))))
    mu.T <- eta.tw+alpha.new*Y
    #tt<-rnorm(1,mu.T,delta.new)
    #log.p.Tb<-dnorm(max(object$yy$Time),mu.T,delta.new,TRUE)
    Y1<-tail(y[id.i],1)
    time1<-eta.tw+alpha.new*Y1
    log.p.Tb<-dnorm(time1,mu.T,delta.new,TRUE)
    log.p.yb + log.p.Tb + log.p.b
  }

