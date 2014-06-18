DSGEVAR.default <- function(dsgedata,lambda=Inf,p=2,ObserveMat,initialvals,partomats,
                            priorform,priorpars,parbounds,parnames=NULL,
                            optimMethod="Nelder-Mead",optimLower=NULL,optimUpper=NULL,optimControl=list(),
                            IRFs=TRUE,irf.periods=20,scalepar=1,keep=50000,burnin=10000){
  #
  cat('Trying to solve the model with your initial values... ')
  dsgemats1t <- partomats(initialvals)
  dsgesolved1t <- SDSGE(dsgemats1t$A,dsgemats1t$B,dsgemats1t$C,dsgemats1t$D,dsgemats1t$F,dsgemats1t$G,dsgemats1t$H,dsgemats1t$J,dsgemats1t$K,dsgemats1t$L,dsgemats1t$M,dsgemats1t$N)
  StateMats1t <- .DSGEstatespace(dsgesolved1t$N,dsgesolved1t$P,dsgesolved1t$Q,dsgesolved1t$R,dsgesolved1t$S)
  cat('Done. \n')
  #
  parametersTrans <- .DSGEParTransform(initialvals,NULL,priorform,parbounds)
  #
  dsgedataret <- dsgedata
  kdata <- .dsgevardata(dsgedata,p,FALSE)
  dsgedata <- kdata$Y
  #
  dsgemode <- NULL
  cat(' \n', sep="")
  cat('Beginning optimization, ', date(),'. \n', sep="")
  cat('Using Optimisation Method: ',optimMethod,'. \n', sep="")
  if(is.finite(lambda)==TRUE){
    if(optimMethod=="Nelder-Mead"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosterior,method="Nelder-Mead",control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else if(optimMethod=="BFGS"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosterior,method="BFGS",control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else if(optimMethod=="CG"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosterior,method="CG",control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else if(optimMethod=="L-BFGS-B"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosterior,method="L-BFGS-B",lower=optimLower,upper=optimUpper,control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else if(optimMethod=="SANN"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosterior,method="SANN",control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else{
      stop("You have entered an unrecognized optimization method.\n",call.=FALSE)
    }
  }else{
    if(optimMethod=="Nelder-Mead"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosteriorInf,method="Nelder-Mead",control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else if(optimMethod=="BFGS"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosteriorInf,method="BFGS",control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else if(optimMethod=="CG"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosteriorInf,method="CG",control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else if(optimMethod=="L-BFGS-B"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosteriorInf,method="L-BFGS-B",lower=optimLower,upper=optimUpper,control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else if(optimMethod=="SANN"){
      dsgemode <- optim(par=parametersTrans,fn=.DSGEVARLogPosteriorInf,method="SANN",control=optimControl,kdata=kdata,lambda=lambda,p=p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
    }else{
      stop("You have entered an unrecognized optimization method.\n",call.=FALSE)
    }
  }
  #
  ConvCode <- dsgemode$convergence
  ConvReport <- 0
  if(ConvCode==0){
    ConvReport <- "successful completion"
  }else if(ConvCode==1){
    ConvReport <- "maximum number of iterations reached"
  }else if(ConvCode==10){
    ConvReport <- "degeneracy of the Nelder-Mead simplex"
  }else if(ConvCode==51){
    ConvReport <- "warning from L-BFGS-B"
  }else if(ConvCode==52){
    ConvReport <- "error from L-BFGS-B"
  }else{}
  #
  cat('Optimization over, ', date(),'. \n', sep="")
  cat(' \n', sep="")
  cat('Optimizer Convergence Code: ',dsgemode$convergence,'; ',ConvReport,'. \n', sep="")
  cat(' \n', sep="")
  cat('Optimizer Iterations: \n', sep="")
  print(dsgemode$counts)
  #
  parametersMode <- .DSGEParTransform(NULL,dsgemode$par,priorform,parbounds)
  #
  parametersModeHessian <- solve(dsgemode$hessian)
  parametersModeHessian <- diag(parametersModeHessian)
  parametersModeHessian <- sqrt(parametersModeHessian)
  #
  ParActualSEs <- .DSGEParTransform(NULL,dsgemode$par,priorform,parbounds) - .DSGEParTransform(NULL,dsgemode$par - parametersModeHessian,priorform,parbounds)
  #
  parametersMode <- matrix(parametersMode,nrow=1)
  ParActualSEs <- matrix(ParActualSEs,nrow=1)
  parametersModeHessian <- matrix(parametersModeHessian,nrow=1)
  #
  ModeTable <- matrix(NA,nrow=length(dsgemode$par),ncol=2)
  ModeTable[,1] <- parametersMode
  ModeTable[,2] <- ParActualSEs
  #
  colnames(ModeTable) <- c("Estimate","SE")
  if(class(parnames)=="character"){
    rownames(ModeTable) <- parnames
  }
  cat(' \n', sep="")
  cat('Parameter Estimates and Standard Errors (SE) at the Posterior Mode: \n', sep="")
  cat(' \n', sep="")
  print(ModeTable)
  #
  rownames(parametersMode) <- "Parameter:"
  rownames(ParActualSEs) <- "Parameter:"
  rownames(parametersModeHessian) <- "Parameter:"
  if(class(parnames)=="character"){
    colnames(parametersMode) <- parnames
    colnames(parametersModeHessian) <- parnames
    colnames(ParActualSEs) <- parnames
  }
  #
  #
  #
  #
  #
  cat(' \n', sep="")
  cat('Trying to Compute DSGE-VAR Prior at Posterior Mode, ', date(),'. \n', sep="")
  dsgeprior <- .DSGEVARPrior(c(parametersMode),dsgedata,kdata$X,p,ObserveMat,partomats,priorform,priorpars,parbounds)
  cat('Done. \n')
  #
  cat(' \n', sep="")
  cat('Beginning DSGE-VAR run, ', date(),'. \n', sep="")
  DSGEMCMCRes <- 0; VARMCMCRes <- 0
  if(is.finite(lambda)==TRUE){
    DSGEMCMCRes <- .DSGEMCMCDraw(dsgemode,scalepar,keep,burnin,kdata,lambda,p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat,partomats,priorform,priorpars,parbounds,parallel=FALSE)
    VARMCMCRes <- .DSGEVARMCMC(DSGEMCMCRes$parameters,c(parametersMode),kdata,lambda,p,ObserveMat,partomats,priorform,priorpars,parbounds)
  }else{
    DSGEMCMCRes <- .DSGEMCMCDrawInf(dsgemode,scalepar,keep,burnin,kdata,lambda,p,YY=kdata$YY,XY=kdata$XY,XX=kdata$XX,ObserveMat,partomats,priorform,priorpars,parbounds,parallel=FALSE)
    VARMCMCRes <- .DSGEVARMCMCInf(DSGEMCMCRes$parameters,c(parametersMode),kdata,lambda,p,ObserveMat,partomats,priorform,priorpars,parbounds)
  }
  cat('MCMC run finished, ', date(),'. \n', sep="")
  #
  cat('Acceptance Rate: ', DSGEMCMCRes$acceptRate,'. \n', sep="")
  #
  if(class(parnames)=="character"){
    colnames(DSGEMCMCRes$parameters) <- parnames
  }
  #
  PostTable <- matrix(NA,nrow=length(dsgemode$par),ncol=4)
  PostTable[,1] <- parametersMode
  PostTable[,2] <- ParActualSEs
  PostTable[,3] <- apply(DSGEMCMCRes$parameters,2,mean)
  PostTable[,4] <- apply(DSGEMCMCRes$parameters,2,sd)
  #
  colnames(PostTable) <- c("Posterior.Mode","SE.Mode","Posterior.Mean","SE.Posterior")
  if(class(parnames)=="character"){
    rownames(PostTable) <- parnames
  }
  cat(' \n', sep="")
  cat('Parameter Estimates and Standard Errors: \n', sep="")
  cat(' \n', sep="")
  print(PostTable)
  #
  IRFDVs <- NULL; A0Mats <- NULL; 
  IRFDs <- NULL
  dsgemats1t <- NULL; dsgesolved1t <- NULL; StateMats1t <- NULL
  if(IRFs == TRUE){
    cat(' \n')
    cat('Computing DSGE IRFs now... ')
    dsgemats1t <- partomats(c(parametersMode))
    dsgesolved1t <- SDSGE(dsgemats1t$A,dsgemats1t$B,dsgemats1t$C,dsgemats1t$D,dsgemats1t$F,dsgemats1t$G,dsgemats1t$H,dsgemats1t$J,dsgemats1t$K,dsgemats1t$L,dsgemats1t$M,dsgemats1t$N)
    StateMats1t <- .DSGEstatespace(dsgesolved1t$N,dsgesolved1t$P,dsgesolved1t$Q,dsgesolved1t$R,dsgesolved1t$S)
    IRFDs <- array(0,dim=c(irf.periods,ncol(StateMats1t$F),nrow(dsgemats1t$N),keep))
    for(i in 1:keep){
      dsgemats <- partomats(DSGEMCMCRes$parameters[i,])
      dsgesolved <- SDSGE(dsgemats$A,dsgemats$B,dsgemats$C,dsgemats$D,dsgemats$F,dsgemats$G,dsgemats$H,dsgemats$J,dsgemats$K,dsgemats$L,dsgemats$M,dsgemats$N)
      iIRF <- IRF(dsgesolved,sqrt(DSGEMCMCRes$parameters[i,(ncol(DSGEMCMCRes$parameters)-nrow(dsgemats1t$N)+1):ncol(DSGEMCMCRes$parameters)]),irf.periods,varnames=NULL,plot=FALSE,save=FALSE)
      IRFDs[,,,i] <- iIRF$IRFs
    }
    cat('Done. \n') 
    #
    cat(' \n')
    cat('Starting DSGE-VAR IRFs, ',date(),'. \n', sep="")
    #
    IRFDVs <- array(NA,dim=c(ncol(dsgedata),ncol(dsgedata),irf.periods*keep))
    #
    A0Mats <- .DSGEVARIRFMatrices(DSGEMCMCRes$parameters,VARMCMCRes$Sigma,p,ObserveMat,partomats,priorform,priorpars,parbounds)
    IRFDVs <- .Call("DSGEVARIRFs", ncol(dsgedata),ncol(dsgedata)*p,0,keep,irf.periods,VARMCMCRes$Phi,A0Mats,IRFDVs, PACKAGE = "BMR", DUP = FALSE)
    IRFDVs <- IRFDVs$ImpStore
    IRFStore <- array(NA,dim=c(ncol(dsgedata),ncol(dsgedata),irf.periods,keep))
    for(i in 1:keep){
      IRFStore[,,,i] <- IRFDVs[,,((i-1)*irf.periods+1):(i*irf.periods)]
    }
    #
    IRFDVs <- 0
    IRFDVs<-apply(IRFStore,c(3,1,2),sort)
    #
    IRFDVs<-aperm(IRFDVs,c(2,3,1,4)); IRFDVs<-aperm(IRFDVs,c(1,2,4,3))
    #
    cat('DSGEVAR IRFs finished, ', date(),'. \n', sep="")
  }
  #
  dsgevarret <- list(Parameters=DSGEMCMCRes$parameters,Beta=VARMCMCRes$Phi,Sigma=VARMCMCRes$Sigma,DSGEIRFs=IRFDs,DSGEVARIRFs=IRFDVs,lambda=lambda,ModeParamTrans=dsgemode$par,ModeHessian=dsgemode$hessian,scalepar=scalepar,AcceptanceRate=DSGEMCMCRes$acceptRate,ObserveMat=ObserveMat,data=dsgedataret,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds)
  class(dsgevarret) <- "DSGEVAR"
  #
  return(dsgevarret)
}

.dsgevardata <- function(mydata,p,constant){
  #
  Tr <- as.numeric(dim(mydata)[1])
  M <- as.numeric(dim(mydata)[2])
  Tp <- Tr - p
  #
  Yraw <- as.matrix(mydata,ncol=M)
  #
  X <- embed(Yraw,p+1); X <- X[,(M+1):ncol(X)]
  if(constant == TRUE){X<-cbind(rep(1,(Tp)),X)}
  #
  K <- as.numeric(dim(X)[2])
  #
  Y <- Yraw[(p+1):nrow(Yraw),]
  #
  YY <- (1/nrow(Y))*t(Y)%*%Y
  XY <- (1/nrow(Y))*t(X)%*%Y
  XX <- (1/nrow(Y))*t(X)%*%X
  #
  return=list(Y=Y,X=X,M=M,K=K,Yraw=Yraw,Tp=Tp,YY=YY,XY=XY,XX=XX)
}

.DSGEVARLogPosterior <- function(dsgeparTrans,kdata,lambda,p,YY,XY,XX,ObserveMat,partomats,priorform,priorpars,parbounds){
  #
  Y <- kdata$Y; X <- kdata$X
  #
  dsgepar <- .DSGEParTransform(NULL,dsgeparTrans,priorform,parbounds)
  #
  logGPR <- .LGPR(ncol(Y),((1+lambda)*nrow(Y))-ncol(Y)*p,(lambda*nrow(Y))-ncol(Y)*p);
  #
  dsgeprior <- .DSGEVARPrior(dsgepar,Y,X,p,ObserveMat,partomats,priorform,priorpars,parbounds)
  GammaYY <- dsgeprior$GammaYY; GammaXY <- dsgeprior$GammaXY; GammaXX <- dsgeprior$GammaXX
  #
  tau <- lambda/(1+lambda);
  GammaBarYY <- (tau*GammaYY)+((1-tau)*YY)
  GammaBarXY <- (tau*t(GammaXY))+((1-tau)*XY)
  GammaBarXX <- (tau*GammaXX)+((1-tau)*XX)
  #
  logLikelihood <- -.Call("DSGEVARLikelihood", logGPR,XX,GammaYY,GammaXY,GammaXX,GammaBarYY,GammaBarXY,GammaBarXX,lambda,nrow(Y),ncol(Y),p, PACKAGE = "BMR", DUP = FALSE)$logLikelihood
  #
  logPosterior <- .DSGEPriors(dsgepar,dsgeparTrans,priorform,priorpars,parbounds,logLikelihood)
  #
  return(logPosterior)
}

.DSGEVARLogPosteriorInf <- function(dsgeparTrans,kdata,lambda,p,YY,XY,XX,ObserveMat,partomats,priorform,priorpars,parbounds){
  #
  Y <- kdata$Y; X <- kdata$X
  #
  dsgepar <- .DSGEParTransform(NULL,dsgeparTrans,priorform,parbounds)
  #
  dsgeprior <- .DSGEVARPrior(dsgepar,Y,X,p,ObserveMat,partomats,priorform,priorpars,parbounds)
  GammaYY <- dsgeprior$GammaYY; GammaXY <- dsgeprior$GammaXY; GammaXX <- dsgeprior$GammaXX
  #
  logLikelihood <- -.Call("DSGEVARLikelihoodInf", YY,XY,XX,GammaYY,GammaXY,GammaXX,nrow(Y),ncol(Y),p, PACKAGE = "BMR", DUP = FALSE)$logLikelihood
  #
  logPosterior <- .DSGEPriors(dsgepar,dsgeparTrans,priorform,priorpars,parbounds,logLikelihood)
  #
  return(logPosterior)
}

.LGPR <- function(n,a,b){
  logGPR = 0
  for(i in 1:n){
    logGPR = logGPR + lgamma((a-i+1)/2) - lgamma((b-i+1)/2)
  }
  return(logGPR)
}

.DSGEVARPrior <- function(parameters,dsgedata,X,p,ObserveMat,partomats,priorform,priorpars,parbounds){
  #
  dsgemats <- partomats(parameters)
  dsgesolved <- SDSGE(dsgemats$A,dsgemats$B,dsgemats$C,dsgemats$D,dsgemats$F,dsgemats$G,dsgemats$H,dsgemats$J,dsgemats$K,dsgemats$L,dsgemats$M,dsgemats$N)
  #
  StateMats <- .DSGEstatespace(dsgesolved$N,dsgesolved$P,dsgesolved$Q,dsgesolved$R,dsgesolved$S)
  #
  SigmaX <- .Call("DSGEVARPriorC", dsgedata,ObserveMat,StateMats$F,StateMats$G,dsgesolved$N,dsgemats$shocks,p,500, PACKAGE = "BMR", DUP = FALSE)$SigmaX
  #
  GammaYY <- SigmaX[,,1]
  #
  #
  #
  GammaXX <- matrix(0,ncol(dsgedata)*p,ncol(dsgedata)*p)
  GammaXY <- matrix(0,ncol(dsgedata),ncol(dsgedata)*p)
  #
  for(i in 1:p){
    RowInd <- (1:ncol(dsgedata))+((i-1)*ncol(dsgedata))
    for(j in 1:p){
      ColInd <- (1:ncol(dsgedata))+((j-1)*ncol(dsgedata))
      if(i==j){
        SelectMat <- SigmaX[,,1]
      }else if(i<j){
        SelectMat <- SigmaX[,,1+j-i]
      }else{
        SelectMat <- t(SigmaX[,,1+i-j])
      }
      GammaXX[RowInd,ColInd] = SelectMat
    }
    GammaXY[,RowInd] = SigmaX[,,1+i]
  }
  #
  return(list(GammaYY=GammaYY,GammaXX=GammaXX,GammaXY=GammaXY))
}

.DSGEMCMCDraw <- function(dsgeopt,scalepar,keep,burnin,kdata,lambda,p,YY,XY,XX,ObserveMat,partomats,priorform,priorpars,parbounds,parallel=FALSE){
  #
  Draws <- matrix(NA,nrow=(keep+burnin+1),ncol=length(dsgeopt$par))
  #
  Draws[1,] <- dsgeopt$par
  if(parallel==TRUE){
    Draws[1,] <- dsgeopt$par + runif(1,-1,1)*c(sqrt(diag(solve(dsgeopt$hessian))))
  }
  #
  PrevLP <- (-1)*dsgeopt$value
  if(parallel==TRUE){
    PrevLP <- (-1)*.DSGEVARLogPosterior(c(Draws[1,]),kdata,lambda,p,YY,XY,XX,ObserveMat,partomats,priorform,priorpars,parbounds)
  }
  #
  PickMeInstead <- matrix(c(Draws[1,]))
  #
  CovM <- scalepar*(solve(dsgeopt$hessian))
  CovMChol <- t(chol(CovM))
  #
  Acceptances <- 0
  #
  for (i in 1:(keep+burnin)){
    proposal <- PickMeInstead + CovMChol%*%matrix(rnorm(length(dsgeopt$par)))
    PropLP <- (-1)*.DSGEVARLogPosterior(c(proposal),kdata,lambda,p,YY,XY,XX,ObserveMat,partomats,priorform,priorpars,parbounds)
    if(is.nan(PropLP)){
      PropLP <- -1000000
    }
    #
    if(runif(1) < exp(PropLP-PrevLP)){
      Draws[i+1,] <- t(proposal)
      if(i > burnin){
        Acceptances <- Acceptances + 1
      }
      PickMeInstead <- proposal
      PrevLP <- PropLP
    }else{
      Draws[i+1,] <- Draws[i,]
    }
  }
  #
  Draws <- Draws[(burnin+2):nrow(Draws),]
  accept <- Acceptances/keep
  #
  for(i in 1:keep){
    Draws[i,] <- .DSGEParTransform(NULL,Draws[i,],priorform,parbounds)
  }
  #
  if(parallel==FALSE){
    return=list(parameters=Draws,acceptRate=accept)
  }else{
    return(list(Draws,accept))
  }
}

.DSGEMCMCDrawInf <- function(dsgeopt,scalepar,keep,burnin,kdata,lambda,p,YY,XY,XX,ObserveMat,partomats,priorform,priorpars,parbounds,parallel=FALSE){
  #
  Draws <- matrix(NA,nrow=(keep+burnin+1),ncol=length(dsgeopt$par))
  #
  Draws[1,] <- dsgeopt$par
  if(parallel==TRUE){
    Draws[1,] <- dsgeopt$par + runif(1,-1,1)*c(sqrt(diag(solve(dsgeopt$hessian))))
  }
  #
  PrevLP <- (-1)*dsgeopt$value
  if(parallel==TRUE){
    PrevLP <- (-1)*.DSGEVARLogPosteriorInf(c(Draws[1,]),kdata,lambda,p,YY,XY,XX,ObserveMat,partomats,priorform,priorpars,parbounds)
  }
  #
  PickMeInstead <- matrix(c(Draws[1,]))
  #
  CovM <- scalepar*(solve(dsgeopt$hessian))
  CovMChol <- t(chol(CovM))
  #
  Acceptances <- 0
  #
  for (i in 1:(keep+burnin)){
    proposal <- PickMeInstead + CovMChol%*%matrix(rnorm(length(dsgeopt$par)))
    PropLP <- (-1)*.DSGEVARLogPosteriorInf(c(proposal),kdata,lambda,p,YY,XY,XX,ObserveMat,partomats,priorform,priorpars,parbounds)
    if(is.nan(PropLP)){
      PropLP <- -1000000
    }
    #
    if(runif(1) < exp(PropLP-PrevLP)){
      Draws[i+1,] <- t(proposal)
      if(i > burnin){
        Acceptances <- Acceptances + 1
      }
      PickMeInstead <- proposal
      PrevLP <- PropLP
    }else{
      Draws[i+1,] <- Draws[i,]
    }
  }
  #
  Draws <- Draws[(burnin+2):nrow(Draws),]
  accept <- Acceptances/keep
  #
  for(i in 1:keep){
    Draws[i,] <- .DSGEParTransform(NULL,Draws[i,],priorform,parbounds)
  }
  #
  if(parallel==FALSE){
    return=list(parameters=Draws,acceptRate=accept)
  }else{
    return(list(Draws,accept))
  }
}

.DSGEVARMCMC <- function(dsgemcmc,dsgemode,kdata,lambda,p,ObserveMat,partomats,priorform,priorpars,parbounds){
  #
  Y <- kdata$Y; X <- kdata$X; YY <- kdata$YY; XY <- kdata$XY; XX <- kdata$XX;
  #
  dsgeprior <- .DSGEVARPrior(dsgemode,Y,X,p,ObserveMat,partomats,priorform,priorpars,parbounds)
  GammaYY <- dsgeprior$GammaYY; GammaXY <- dsgeprior$GammaXY; GammaXX <- dsgeprior$GammaXX
  #
  tau <- lambda/(1+lambda)
  GammaBarYY <- (tau*GammaYY)+((1-tau)*YY)
  GammaBarXY <- (tau*t(GammaXY))+((1-tau)*XY)
  GammaBarXX <- (tau*GammaXX)+((1-tau)*XX)
  #
  GammaBarYY <- array(0,dim=c(c(dim(YY)),nrow(dsgemcmc)))
  GammaBarXY <- array(0,dim=c(c(dim(XY)),nrow(dsgemcmc)))
  GammaBarXX <- array(0,dim=c(c(dim(XX)),nrow(dsgemcmc)))
  GXX <- array(0,dim=c(c(dim(XX)),nrow(dsgemcmc)))
  #
  for(t in 1:nrow(dsgemcmc)){
    dsgeparameters <- dsgemcmc[t,]
    dsgeprior <- .DSGEVARPrior(dsgeparameters,Y,X,p,ObserveMat,partomats,priorform,priorpars,parbounds)
    GammaYY <- dsgeprior$GammaYY; GammaXY <- dsgeprior$GammaXY; GammaXX <- dsgeprior$GammaXX
    #
    GammaBarYY[,,t] <- (tau*GammaYY)+((1-tau)*YY)
    GammaBarXY[,,t] <- (tau*t(GammaXY))+((1-tau)*XY)
    GammaBarXX[,,t] <- (tau*GammaXX)+((1-tau)*XX)
    #
    GXX[,,t] <- GammaXX
    #
  }
  #
  RepsRun <- .Call("DSGEVARReps", GammaBarYY,GammaBarXY,GammaBarXX,GXX,XX,lambda,nrow(dsgemcmc),nrow(Y),ncol(Y),p, PACKAGE = "BMR", DUP = FALSE)
  #
  return(list(Phi=RepsRun$Beta,Sigma=RepsRun$Sigma))
}

.DSGEVARMCMCInf <- function(dsgemcmc,dsgemode,kdata,lambda,p,ObserveMat,partomats,priorform,priorpars,parbounds){
  #
  Y <- kdata$Y; X <- kdata$X
  #
  dsgeprior <- .DSGEVARPrior(dsgemode,Y,X,p,ObserveMat,partomats,priorform,priorpars,parbounds)
  GammaYY <- dsgeprior$GammaYY; GammaXY <- dsgeprior$GammaXY; GammaXX <- dsgeprior$GammaXX
  #
  GammaBarYY <- array(0,dim=c(c(dim(GammaYY)),nrow(dsgemcmc)))
  GammaBarXY <- array(0,dim=c(c(dim(t(GammaXY))),nrow(dsgemcmc)))
  GammaBarXX <- array(0,dim=c(c(dim(GammaXX)),nrow(dsgemcmc)))
  #
  for(t in 1:nrow(dsgemcmc)){
    dsgeparameters <- dsgemcmc[t,]
    dsgeprior <- .DSGEVARPrior(dsgeparameters,Y,X,p,ObserveMat,partomats,priorform,priorpars,parbounds)
    GammaYY <- dsgeprior$GammaYY; GammaXY <- dsgeprior$GammaXY; GammaXX <- dsgeprior$GammaXX
    #
    GammaBarYY[,,t] <- GammaYY
    GammaBarXY[,,t] <- t(GammaXY)
    GammaBarXX[,,t] <- GammaXX
    #
  }
  #
  RepsRun <- .Call("DSGEVARRepsInf", GammaBarYY,GammaBarXY,GammaBarXX,lambda,nrow(dsgemcmc),nrow(Y),ncol(Y),p, PACKAGE = "BMR", DUP = FALSE)
  #
  return(list(Phi=RepsRun$Beta,Sigma=RepsRun$Sigma))
}

.DSGEVARIRFMatrices <- function(dsgepars,Sigma,p,ObserveMat,partomats,priorform,priorpars,parbounds){
  #
  A0Mats <- array(0,dim=dim(Sigma))
  #
  parameters <- c(dsgepars[1,]); A0 <- matrix(0,dim(Sigma)[1],dim(Sigma)[2])
  #
  for(i in 1:(dim(Sigma)[3])){
    parameters <- c(dsgepars[i,])
    SigmaEpsilon <- Sigma[,,i]
    #
    dsgemats <- partomats(parameters)
    dsgesolved <- SDSGE(dsgemats$A,dsgemats$B,dsgemats$C,dsgemats$D,dsgemats$F,dsgemats$G,dsgemats$H,dsgemats$J,dsgemats$K,dsgemats$L,dsgemats$M,dsgemats$N)
    #
    StateMats <- .DSGEstatespace(dsgesolved$N,dsgesolved$P,dsgesolved$Q,dsgesolved$R,dsgesolved$S)
    #
    Shocks <- matrix(0,nrow(StateMats$G),nrow(StateMats$G))
    Shocks[(nrow(dsgemats$shocks)+1):nrow(StateMats$G),(nrow(dsgemats$shocks)+1):nrow(StateMats$G)] <- sqrt(dsgemats$shocks)
    #
    Shocks <- StateMats$G%*%Shocks
    #
    GMatShocks <- Shocks[,(nrow(StateMats$G)-nrow(dsgemats$shocks)+1):nrow(StateMats$G)]
    #
    SigmaChol <- t(chol(SigmaEpsilon))
    #
    AMat = t(ObserveMat)%*%GMatShocks
    #
    QRAMat <- qr(t(AMat))
    Q <- (-1)*qr.Q(QRAMat)
    R <- (-1)*qr.R(QRAMat)
    #
    RVec <- diag(R)
    PosElem <- (RVec>0)*matrix(1:ncol(ObserveMat))
    PosElem <- PosElem[PosElem>0,]
    S <- matrix(-1,ncol(ObserveMat),1)
    if(length(PosElem)>0){
      S[PosElem,] <- 1
    }
    S <- diag(c(S))
    Q = Q%*%S
    A0 = SigmaChol%*%t(Q)
    #
    A0Mats[,,i] <- A0
    #
  }
  #
  return(A0Mats)
}