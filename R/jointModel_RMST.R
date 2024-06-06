#' Title Joint model for Longitudinal and Survival Data based on restricted mean survival time (RMST)
#' @description
#' This function fits shared parameter models for the joint modelling of normal longitudinal responses
#' and time-to-event data
#' @usage jointModelRMST(lmeObjects,survObjects,Time,timeVar)
#' @param lmeObjects an object inheriting from class lme
#' @param survObjects an object inheriting from class geeglm
#' @param Time Vector representing the survival time in the data
#' @param timeVar a character string indicating the time variable in the linear mixed effects model
#'
#' @return A list of class jointModel_RMST with components:
#' \describe{
#'  \item{coefficients} {a list with the estimated coefficients. The components of this list are:
#' \item{betas} {the vector of fixed effects for the linear mixed effects model.}
#' \item{sigma}{the measurement error standard deviation for the linear mixed effectsmodel.}
#' \item{gammas}{the vector of baseline covariates for the survival model.}
#' \item{alpha}{the association parameter(s).}
#' \item{delta}{the measurement error standard deviation for the survival model.}
#' \item{D}{the variance-covariance matrix of the random effects.}}
#' \item{Hessian}{the Hessian matrix evaluated at the estimated parameter values.}
#' \item{logLik}{the log-likelihood value.}
#' \item{EB}{a list with components:
#' \item{post.b}{the estimated random effects values.}
#' \item{post.vb}{the estimated variance for the random effects estimates.}
#' \item{Zb}{the estimated random effects part of the linear predictor for the longitudinal outcome}
#' \item{Ztimeb}{the estimated random effects part of the linear predictor for the survival outcome}}
#' \item{n}{the number of subjects}
#' \item{N}{the total number of repeated measurements for the longitudinal outcome.}
#' \item{ni}{a vector with the number of repeated measurements for each subject.}
#' \item{id}{the grouping vector for the longitudinal responses.}
#' \item{x}{a list with the design matrices for the longitudinal and survival model.}
#' \item{yy}{a list with the response vectors for the longitudinal and survival model.}
#' \item{times}{the value of the timeVar argument}
#' \item{data}{a data.frame containing the whole datasetes}
#' \item{data.id}{a data.frame containing the variables for the linear mixed effects model at the time of the event.}
#' \item{termsYx}{the terms component  for the fixed effects of the lmeObjects.}
#' \item{termsYz}{the terms component  for the random effects of the lmeObjects.}
#' \item{termsT}{the terms component of the survObjects.}
#' \item{formYx}{the formula for the fixed effects part of the longitudinal model.}
#' \item{formYz}{the formula for the random effects part of the longitudinal model.}
#' \item{formT}{the formula for the survival model.}
#' \item{timeVar}{a character string indicating the time variable in the linear mixed effects model}
#' \item{call}{the matched call.}}
#' @export
#' @references Zhang C, Li Z, Yang Z, Huang B, Hou Y, Chen Z*. A dynamic prediction model supporting individual life expectancy prediction based on longitudinal time-dependent covariates. IEEE Journal of Biomedical and Health Informatics. 2023. 27(9): 4623-4632
#'
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
#'result<-data.frame(beta=round(c(summary(object5)[[1]][,1],summary(object5)[[2]][,1]),3),sd=round(c(summary(object5)[[1]][,2],summary(object5)[[2]][,2]),3),
#'                   low=round(c(summary(object5)[[1]][,4],summary(object5)[[2]][,4]),3),upp=round(c(summary(object5)[[1]][,5],summary(object5)[[2]][,5]),3),
#'                   p=round(c(summary(object5)[[1]][,6],summary(object5)[[2]][,6]),3))
jointModelRMST<-function(lmeObjects,survObjects,Time,timeVar){
  #initial values
  #####longdinal
  cl <- match.call()
  idT<-lmeObjects$groups[[1]]
  id<-match(idT,unique(idT))
  b<-data.matrix(ranef(lmeObjects))
  D<-lapply(pdMatrix(lmeObjects$modelStruct$reStruct),"*",lmeObjects$sigma^2)[[1]]
  TermsX<-lmeObjects$terms
  data <- lmeObjects$data[all.vars(TermsX)]#
  data <- data[complete.cases(data), ]
  mfX<-model.frame(TermsX,data=data)
  y<-model.response(mfX,"numeric")
  formYx<-formula(lmeObjects)
  X<-model.matrix(formYx,mfX)#固定效应
  formYz<-formula(lmeObjects$modelStruct$reStruct[[1]])
  mfz<-model.frame(terms(formYz),data=data)
  TermsZ<-attr(mfz,"terms")
  Z<-model.matrix(formYz,mfz)#随机效应
  data.id <- data[!duplicated(id), ]#
  data.id[[timeVar]]<-Time
  Xtime<-model.matrix(formula(lmeObjects),model.frame(TermsX,data=data.id))#固定效应
  Ztime<-model.matrix(formYz,model.frame(terms(formYz),data=data.id))#随机效应
  long<-c(X %*% fixef(lmeObjects))+rowSums(Z*b[id,])
  long.id<-tapply(long, id, tail,1)
  ######Pseudo RMST
  W<-survObjects$geese$X
  WW<-W[,-1]
  p1<-survObjects$y
  n<-length(p1)
  N<-length(y)
  ni <- as.vector(tapply(id, id, length))
  fit<-geeglm(p1~WW+long.id,id=c(1:length(p1)))
  #intercept<-fit$coefficients["(Intercept)"]
  #gamma<-fit$coefficients[1:ncol(WW)+1]
  #gammas<-c(intercept,gamma)
  gammas<-head(fit$coefficients,-1)
  alpha<-fit$coefficients["long.id"]
  betas<-fixef(lmeObjects)
  sigma<-lmeObjects$sigma
  #delta<-sqrt(var(fit$residuals))
  delta<-sqrt(sum(fit$residuals^2)/(n-(ncol(W)+1)))
  initial_values<-c(list(betas=betas,sigma=sigma,D=D,gammas=gammas,alpha=alpha,delta=delta))
  x<-list(X=X,Z=Z,W=W,idT=idT,Xtime=Xtime,Ztime=Ztime)
  yy<-list(y=y,p1=p1,Time=Time)
  XtX <- crossprod(X)
  ZtZ <- lapply(split(Z, id), function (x) crossprod(matrix(x, ncol = ncol(Z))))
  names(ZtZ) <- NULL
  ZtZ <- matrix(unlist(ZtZ), n, ncol(Z) * ncol(Z), TRUE)
  ncz=ncol(Z)
  ##Gauss-Hermite quadrature rule components
  GHk<-if(ncz<3 && nrow(Z)<2000) 15 else 9
  GH <- JM:::gauher(GHk)#11个点
  b <- as.matrix(expand.grid(rep(list(GH$x), ncol(Z))))
  k <- nrow(b)
  wGH <- as.matrix(expand.grid(rep(list(GH$w), ncol(Z))))
  wGH <- 2^(ncol(Z)/2) * apply(wGH, 1, prod) * exp(rowSums(b * b))
  diag.D <- !is.matrix(D)
  if (!diag.D) dimnames(D) <- NULL else names(D) <- NULL
  inv_chol_D<-solve(chol(solve(D)))
  det_inv_chol_D<-det(inv_chol_D)
  b <- sqrt(2) * t(inv_chol_D %*% t(b))
  wGH <- wGH * det_inv_chol_D
  dimnames(b) <- NULL
  b2 <- if (ncol(Z) == 1) b * b else t(apply(b, 1, function (x) x %o% x))
  Ztb <- Z %*% t(b)
  Ztime.b <- Ztime %*% t(b)
  environment(opt.survRMST) <- environment(gr.survRMST) <- environment(gr.longRMST) <- environment()
  environment(LogLik.RMSTGH) <- environment(Score.RMSTGH) <- environment()
  ##EM iterations
  iter=50
  Y.mat<-matrix(0,iter+1,ncol(X)+1)
  T.mat<-matrix(0,iter+1,ncol(W)+2)
  B.mat<-matrix(0,iter+1,ncol(Z)*ncol(Z))
  lgLik<-numeric(iter+1)
  conv<-TRUE
  verbose=FALSE
  parscale=NULL
  optimizer="optim"
  iter.qN=350
  only.EM=FALSE
  for(it in 1:iter){
    #save parameter values in matrix
    Y.mat[it, ] <- c(betas, sigma)
    T.mat[it, ] <- c(gammas, alpha,delta)
    B.mat[it,] <- D
    # linear predictors
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    #E-step
    mu.y<-eta.yx+Ztb
    logNorm<-dnorm(y,mu.y,sigma,TRUE)
    log.p.yb<-rowsum(logNorm,id,reorder=FALSE)
    dimnames(log.p.yb)<-NULL
    mu.T<-eta.t
    log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
    dimnames(log.p.Tb)<-NULL
    log.p.b<-rep(JM:::dmvnorm(b,rep(0,ncol(Z)),D,TRUE),each=n)
    p.ytb<-exp(log.p.yb+log.p.Tb+log.p.b)
    p.yt<-c(p.ytb %*% wGH)
    p.byt<-p.ytb/p.yt
    post.b<-p.byt %*% (b*wGH)
    post.vb<-if (ncol(Z)==1){
      c(p.byt %*% (b2*wGH))-c(post.b*post.b)
    }else{
      p.byt %*%(b2*wGH)-t(apply(post.b,1,function(x) x %o% x))
    }
    ##compute log-likelihood
    log.p.yt<-log(p.yt)
    lgLik[it]<-sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
    # print results if verbose
    if (verbose) {
      cat("\n\niter:", it, "\n")
      cat("log-likelihood:", lgLik[it], "\n")
      cat("betas:", round(betas, 4), "\n")
      cat("sigma:", round(sigma, 4), "\n")
      if (!is.null(WW))
        cat("gammas:", round(gammas, 4), "\n")
      cat("alpha:", round(alpha, 4), "\n")
      cat("D:", if (!diag.D) round(D[lower.tri(D, TRUE)], 4) else round(D, 4), "\n")
    }
    # check convergence
    if (it > 5) {
      if (lgLik[it] > lgLik[it - 1]) {
        thets1 <- c(Y.mat[it - 1, ], T.mat[it - 1, ], B.mat[it - 1, ])
        thets2 <- c(Y.mat[it, ], T.mat[it, ], B.mat[it, ])
        check1 <- max(abs(thets2 - thets1) / (abs(thets1) + 1e-03)) < 1e-04
        check2 <- (lgLik[it] - lgLik[it - 1]) < 1e-10 * (abs(lgLik[it - 1]) + 1e-10)
        if (check1 || check2) {
          conv <- FALSE
          if (verbose)
            cat("\n\nconverged!\n")
          break}
      }
    }
    if (iter == 0) break
    #M-step
    Zb<-rowSums(Z*post.b[id,],na.rm=TRUE)
    mu<-y-eta.yx
    tr.tZZvarb<-sum(ZtZ*post.vb,na.rm=TRUE)
    sigman<-sqrt(c(crossprod(mu,mu-2*Zb)+crossprod(Zb)+tr.tZZvarb)/N)
    Dn<-matrix(colMeans(p.byt %*% (b2*wGH),na.rm=TRUE),ncol(Z),ncol(Z))
    Dn<-if (diag.D) diag(Dn) else 0.5 * (Dn + t(Dn))
    Hbetas<-JM:::nearPD(JM:::fd.vec(betas,gr.longRMST))
    scbetas<-gr.longRMST(betas)
    betasn<-betas-c(solve(Hbetas,scbetas))

    list.thetas <- list(gammas = gammas, alpha = alpha,log.delta=log(delta))
    list.thetas <- list.thetas[!sapply(list.thetas, is.null)]
    thetas <- unlist(as.relistable(list.thetas))
    optz.surv <- optim(thetas, opt.survRMST, gr.survRMST, method = "BFGS",
                       control = list(maxit = if (it < 5) 20 else 5,
                                      parscale = if (it < 5) rep(0.01, length(thetas)) else rep(0.1, length(thetas))))
    thetasn <- relist(optz.surv$par, skeleton = list.thetas)

    # update parameter values

    betas <- betasn
    sigma <- sigman
    D <- Dn
    gammas <- thetasn$gammas
    alpha <- thetasn$alpha
    delta<-exp(thetasn$log.delta)

  }
  list.thetas <- list(betas = betas, log.sigma = log(sigma), gammas = gammas, alpha = alpha, log.delta=log(delta), D = if (diag.D) log(D) else JM:::chol.transf(D))

  list.thetas <- list.thetas[!sapply(list.thetas, is.null)]

  thetas <- unlist(as.relistable(list.thetas))
  lgLik <- - LogLik.RMSTGH(thetas)
  # if not converged, start quasi-Newton iterations
  if (conv && !only.EM) {
    if (is.null(parscale))
      parscale <- rep(0.01, length(thetas))
    if (verbose)
      cat("\n\nquasi-Newton iterations start.\n\n")
    out <- if (optimizer == "optim") {
      optim(thetas, LogLik.RMSTGH, Score.RMSTGH, method = "BFGS",
            control = list(maxit =iter.qN, parscale = parscale,
                           trace = 10 * verbose))
    } else {
      nlminb(thetas, LogLik.RMSTGH,Score.RMSTGH, scale =parscale,
             control = list(iter.max = iter.qN, trace = 1 * verbose))
    }
    if ((conv <- out$convergence) == 0 || - out[[2]] > lgLik) {
      lgLik <- - out[[2]]
      thetas <- relist(out$par, skeleton = list.thetas)
      betas <- thetas$betas
      sigma <- exp(thetas$log.sigma)
      gammas <- thetas$gammas
      alpha <- thetas$alpha
      delta<- exp(thetas$log.delta)
      D <- thetas$D
      D <- if (diag.D) exp(D) else JM:::chol.transf(D)
      it <- it + if (optimizer == "optim") out$counts[1] else out$iterations
      # compute posterior moments for thetas after quasi-Newton
      eta.yx <- as.vector(X %*% betas)
      eta.tw <- as.vector(W %*% gammas)
      Y <- as.vector(Xtime %*% betas) + Ztime.b
      eta.t <- eta.tw + alpha * Y
      mu.y <- eta.yx + Ztb
      logNorm <- dnorm(y, mu.y, sigma, TRUE)
      log.p.yb <- rowsum(logNorm, id)
      mu.T<-eta.t
      log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
      log.p.b<-rep(JM:::dmvnorm(b,rep(0,ncol(Z)),D,TRUE),each=n)
      p.ytb <- exp(log.p.yb + log.p.Tb + log.p.b)
      p.yt <- c(p.ytb %*% wGH)
      p.byt <- p.ytb / p.yt
      post.b <- p.byt %*% (b * wGH)
      post.vb <- if (ncol(Z) == 1) {
        c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
      } else {
        (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
      }
      Zb <- if (ncol(Z) == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
    }
  }
  Score <- Score.RMSTGH(unlist(thetas))
  # calculate Hessian matrix
  if (verbose) cat("\ncalculating Hessian...\n")
  numeriDeriv = "fd"
  Hessian <- if (numeriDeriv == "fd") {
    JM:::fd.vec(unlist(thetas), Score.RMSTGH, eps = 1e-06)
  } else {
    JM:::cd.vec(unlist(thetas), Score.RMSTGH, eps = 1e-06)
  }
  names(betas) <- names(initial_values$betas)
  if (!diag.D) dimnames(D) <- dimnames(initial_values$D) else names(D) <- names(initial_values$D)
  names(gammas) <- colnames(W)
  nams <- c(paste("Y.", c(names(betas), "sigma"), sep = ""), paste("T.", c(names(gammas), "alpha","delta"), sep = ""),
            paste("B.", if (!diag.D) paste("D", seq(1, ncz * (ncz + 1) / 2), sep = "") else names(D), sep = ""))
  dimnames(Hessian) <- list(nams, nams)
  colnames(post.b) <- colnames(Z)
  out<-list(coefficients = list(betas = betas, sigma = sigma, gammas = gammas, alpha = alpha,
                                delta=delta, D = as.matrix(D)), Hessian = Hessian,
            logLik = lgLik, EB = list(post.b = post.b, post.vb = post.vb,
                                      Zb = if (iter == 0) rowSums(Z * post.b[id, ], na.rm = TRUE) else Zb,
                                      Ztimeb = rowSums(Ztime * post.b)),  iters = it, convergence = conv,
            n = n, N = N, ni = ni, id = id)
  H <- out$Hessian
  if (any(is.na(H) | !is.finite(H))) {
    warning("infinite or missing values in Hessian at convergence.\n")
  } else {
    ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
    if (!all(ev >= -1e-06 * abs(ev[1])))
      warning("Hessian matrix at convergence is not positive definite.\n")
  }
  out$coefficients <- out$coefficients[!sapply(out$coefficients, is.null)]
  out$x <- x
  out$yy <- yy
  out$times <- data[[timeVar]]
  out$data <- data
  out$data.id <- data.id
  out$termsYx <- TermsX
  out$termsYz <- TermsZ
  out$termsT <- survObjects$terms
  out$formYx <- formYx
  out$formYz <- formYz
  out$formT <- survObjects$geese$formula
  out$timeVar <- timeVar
  #out$assignY <- attr(lmeObjects$fixDF, "assign")[-1]
  #out$assignT <- survObjects$assign
  out$call <- cl
  class(out) <- "jointModel_RMST"
  out}


gr.longRMST <-
  function (betas) {
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + alpha * Y
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    sc2<-numeric(ncol(X))
    for (i in 1:ncol(X)){
      #sc2[i]<--sum(alpha*Xtime[,i]/delta^2*c((p.byt*(p1-eta.t)) %*% wGH ),na.rm=TRUE)
      sc2[i]<--sum(alpha*Xtime[,i]/delta^2*p1-alpha*Xtime[,i]/delta^2*c((p.byt*eta.t)  %*% wGH),na.rm=TRUE)
    }
    c(sc1 + sc2)
  }
gr.survRMST <-
  function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    delta<-exp(thetas$log.delta)
    #gammas <- thetas[seq_len(ncol(W))]
    #alpha <- thetas[(ncol(W)+1)]
    #delta<-thetas[(ncol(W)+2)]
    eta.yxT <- as.vector(Xtime %*% betas)
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
    Y <- eta.yxT + Ztime.b
    eta.t <- eta.tw + alpha * Y
    sc.gammas<-numeric(ncol(W))
    for (i in 1:ncol(W))
    {
      sc.gammas[i]<- - sum(W[,i]/delta^2 * p1 -  W[,i]/delta^2*c((p.byt*eta.t)%*% wGH), na.rm = TRUE)
    }
    sc.alpha <- - sum((p.byt*(p1-eta.t)*Y/delta^2)%*% wGH, na.rm = TRUE)
    #sc.alpha <- - sum((p.byt*Y/delta^2 * (p1 -  eta.t))%*% wGH, na.rm = TRUE)
    sc.delta<--delta*(-n/delta+sum((p.byt *((p1-eta.t)^2/delta^3))%*% wGH,na.rm = TRUE))
    #sc.delta<--(-n/delta+sum((p.byt *((p1-eta.t)^2/delta^3))%*% wGH,na.rm = TRUE))
    c(sc.gammas, sc.alpha,sc.delta)
  }
opt.survRMST <-
  function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    delta<-exp(thetas$log.delta)
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else 0
    eta.t <-  eta.tw + alpha * Y
    mu.T<-eta.t
    log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
    p.bytn <- p.byt * log.p.Tb
    -sum(p.bytn %*% wGH, na.rm = TRUE)
  }
Score.RMSTGH <-
  function (thetas) {
    #betas <- thetas[1:ncol(X)]
    #sigma <- exp(thetas[ncol(X) + 1])
    #gammas <- thetas[seq(ncol(X) + 2, ncol(X) + 1 + ncol(W))]
    #alpha <- thetas[ncol(X) + ncol(W) + 2]
    #delta<-exp(thetas[ncol(X)+ncol(W)+3])
    #D <- thetas[seq(ncol(X) + ncol(W) +4, length(thetas))]
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    delta<-exp(thetas$log.delta)
    D <- thetas$D
    D <- if (diag.D) exp(D) else JM:::chol.transf(D)
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id); dimnames(log.p.yb) <- NULL
    mu.T<-eta.t
    log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
    log.p.b<-rep(JM:::dmvnorm(b,rep(0,ncol(Z)),D,TRUE),each=n)
    p.ytb <- exp(log.p.yb + log.p.Tb + log.p.b)
    dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    p.byt <- p.ytb / p.yt
    post.b <-  p.byt %*% (b * wGH)
    post.vb <- if (ncol(Z) == 1) {
      c(p.byt %*% (b2 * wGH)) - c(post.b * post.b)
    } else {
      (p.byt %*% (b2 * wGH)) - t(apply(post.b, 1, function (x) x %o% x))
    }
    Zb <- if (ncol(Z) == 1) post.b[id] else rowSums(Z * post.b[id, ], na.rm = TRUE)
    mu <- y - eta.yx
    tr.tZZvarb <- sum(ZtZ * post.vb, na.rm = TRUE)
    sc1 <- - crossprod(X, y - eta.yx - Zb) / sigma^2
    sc2<-numeric(ncol(X))
    for (i in 1:ncol(X)){
      sc2[i]<--sum(alpha*Xtime[,i]/delta^2*p1-alpha*Xtime[,i]/delta^2*c((p.byt*eta.t)  %*% wGH),na.rm=TRUE)
    }
    score.y <- c(c(sc1 + sc2),
                 - sigma * (- N / sigma + drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    #score.y <- c(c(sc1 + sc2),
    #- 1/(2*sigma) * (- N / sigma + drop(crossprod(mu, mu - 2 * Zb) + crossprod(Zb) + tr.tZZvarb) / sigma^3))
    sc.gammas<-numeric(ncol(W))
    for (i in 1:ncol(W))
    {
      sc.gammas[i]<- - sum(W[,i]/delta^2 * p1 -  W[,i]/delta^2*c((p.byt*eta.t)%*% wGH), na.rm = TRUE)
    }
    sc.alpha <- - sum((p.byt*(p1-eta.tw-alpha*Y)*Y/delta^2)%*% wGH, na.rm = TRUE)
    sc.delta<--delta*(-n/delta+sum((p.byt *((p1-eta.t)^2/delta^3))%*% wGH,na.rm = TRUE))
    #sc.delta<--(-n/delta+sum((p.byt *((p1-eta.t)^2/delta^3))%*% wGH,na.rm = TRUE))
    score.t <- c(sc.gammas, sc.alpha,sc.delta)
    score.b <- if (diag.D) {
      svD <- 1 / D
      svD2 <- svD^2
      cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
      dim(cS.postVB) <- c(ncz, ncz)
      D * 0.5 * (n * svD - diag(cS.postVB) * svD2 - colSums(as.matrix(post.b^2), na.rm = TRUE) * svD2)
    } else {
      svD <- solve(D)
      dD <- JM:::deriv.D(D)
      ndD <- length(dD)
      D1 <- sapply(dD, function (x) sum(svD * x))
      D2 <- t(sapply(dD, function (x) c(svD %*% x %*% svD)))
      cS.postVB <- colSums(as.matrix(post.vb), na.rm = TRUE)
      out <- numeric(ndD)
      for (i in seq_along(dD)) {
        D.mat <- D2[i, ]
        dim(D.mat) <- c(ncz, ncz)
        out[i] <- sum(D2[i, ] * cS.postVB, na.rm = TRUE) + sum((post.b %*% D.mat) * post.b, na.rm = TRUE)
      }
      J <- JM:::jacobian2(attr(D, "L"), ncz)
      drop(0.5 * (n * D1 - out) %*% J)
    }
    c(score.y, score.t, score.b)
  }
LogLik.RMSTGH <-
  function (thetas) {
    thetas <- relist(thetas, skeleton = list.thetas)
    betas <- thetas$betas
    sigma <- exp(thetas$log.sigma)
    gammas <- thetas$gammas
    alpha <- thetas$alpha
    delta<-exp(thetas$log.delta)
    D <- thetas$D
    D <- if (diag.D) exp(D) else JM:::chol.transf(D)
    # linear predictors
    eta.yx <- as.vector(X %*% betas)
    eta.yxT <- as.vector(Xtime %*% betas)
    Y <- eta.yxT + Ztime.b
    eta.tw <- if (!is.null(W)) as.vector(W %*% gammas) else rep(0, n)
    eta.t <- eta.tw + alpha * Y
    mu.y <- eta.yx + Ztb
    logNorm <- dnorm(y, mu.y, sigma, TRUE)
    log.p.yb <- rowsum(logNorm, id); dimnames(log.p.yb) <- NULL
    mu.T<-eta.t
    log.p.Tb<-dnorm(p1,mu.T,delta,TRUE)
    log.p.b<-rep(JM:::dmvnorm(b,rep(0,ncol(Z)),D,TRUE),each=n)
    p.ytb <- exp(log.p.yb + log.p.Tb + log.p.b)
    dimnames(p.ytb) <- NULL
    p.yt <- c(p.ytb %*% wGH)
    p.byt <- p.ytb / p.yt
    log.p.yt <- log(p.yt)
    - sum(log.p.yt[is.finite(log.p.yt)], na.rm = TRUE)
  }



