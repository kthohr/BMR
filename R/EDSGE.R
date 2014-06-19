EDSGE.default <- function(dsgedata,chains=1,cores=1,ObserveMat,initialvals,partomats,
                          priorform,priorpars,parbounds,parnames=NULL,
                          optimMethod="Nelder-Mead",optimLower=NULL,optimUpper=NULL,optimControl=list(),
                          DSGEIRFs=TRUE,irf.periods=20,scalepar=1,keep=50000,burnin=10000){
  #
  cat('Trying to solve the model with your initial values... ')
  dsgemats1t <- partomats(initialvals)
  dsgesolved1t <- SDSGE(dsgemats1t$A,dsgemats1t$B,dsgemats1t$C,dsgemats1t$D,dsgemats1t$F,dsgemats1t$G,dsgemats1t$H,dsgemats1t$J,dsgemats1t$K,dsgemats1t$L,dsgemats1t$M,dsgemats1t$N)
  StateMats1t <- .DSGEstatespace(dsgesolved1t$N,dsgesolved1t$P,dsgesolved1t$Q,dsgesolved1t$R,dsgesolved1t$S)
  cat('Done. \n')
  #
  parametersTrans <- .DSGEParTransform(initialvals,NULL,priorform,parbounds)
  #
  dsgemode <- NULL
  cat(' \n', sep="")
  cat('Beginning optimization, ', date(),'. \n', sep="")
  cat('Using Optimisation Method: ',optimMethod,'. \n', sep="")
  if(optimMethod=="Nelder-Mead"){
    dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="Nelder-Mead",control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
  }else if(optimMethod=="BFGS"){
    dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="BFGS",control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
  }else if(optimMethod=="CG"){
    dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="CG",control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
  }else if(optimMethod=="L-BFGS-B"){
    dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="L-BFGS-B",lower=optimLower,upper=optimUpper,control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
  }else if(optimMethod=="SANN"){
    dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="SANN",control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
  }else{
    stop("You have entered an unrecognized optimization method.\n",call.=FALSE)
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
  logMargLikelihood <- .LaplaceMargLikelihood(dsgemode)
  #
  cat(' \n', sep="")
  cat('Log Marginal Likelihood: ',logMargLikelihood,'. \n', sep="")
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
  cat(' \n', sep="")
  cat('Beginning MCMC run, ', date(),'. \n', sep="")
  MCMCRes <- 0
  if(chains==1){
    MCMCRes <- .DSGEMCMC(dsgemode,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
  }else if(chains > 1){
    MCMCRes <- .DSGEMCMCMulti(dsgemode,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds,chains,cores)
  }
  cat('MCMC run finished, ', date(),'. \n', sep="")
  if(chains==1){
    cat('Acceptance Rate: ', MCMCRes$acceptRate,'. \n', sep="")
  }else{
    cat('Acceptance Rate: ', sep="")
    for(kk in 1:(chains-1)){
      cat('Chain ',kk,': ', MCMCRes$acceptRate[kk],'; ', sep="")
    }
    cat('Chain ',chains,': ', MCMCRes$acceptRate[chains],'. \n', sep="")
    #
    # Chain convergence statistics:
    #
    Diagnostics <- matrix(.MCMCDiagnostics(MCMCRes$parameters,chains),nrow=1)
    #
    rownames(Diagnostics) <- "Stat:"
    if(class(parnames)=="character"){
      colnames(Diagnostics) <- parnames
    }
    #
    cat(' \n', sep="")
    cat('Root-R Chain-Convergence Statistics: \n', sep="")
    print(Diagnostics)
    cat(' \n', sep="")
  }
  #
  if(class(parnames)=="character"){
    colnames(MCMCRes$parameters) <- parnames
  }
  #
  PostTable <- matrix(NA,nrow=length(dsgemode$par),ncol=4)
  PostTable[,1] <- parametersMode
  PostTable[,2] <- ParActualSEs
  PostTable[,3] <- apply(MCMCRes$parameters,2,mean)
  PostTable[,4] <- apply(MCMCRes$parameters,2,sd)
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
  IRFs <- NULL
  if(DSGEIRFs == TRUE){
    cat(' \n')
    cat('Computing IRFs now... ')
    IRFs <- array(0,dim=c(irf.periods,ncol(StateMats1t$F),nrow(dsgemats1t$N),keep))
    for(i in 1:keep){
      dsgemats <- partomats(MCMCRes$parameters[i,])
      dsgesolved <- SDSGE(dsgemats$A,dsgemats$B,dsgemats$C,dsgemats$D,dsgemats$F,dsgemats$G,dsgemats$H,dsgemats$J,dsgemats$K,dsgemats$L,dsgemats$M,dsgemats$N)
      iIRF <- IRF(dsgesolved,sqrt(MCMCRes$parameters[i,(ncol(MCMCRes$parameters)-nrow(dsgemats1t$N)+1):ncol(MCMCRes$parameters)]),irf.periods,varnames=NULL,plot=FALSE,save=FALSE)
      IRFs[,,,i] <- iIRF$IRFs
    }
    cat('Done. \n')
  }
  #
  dsgeret <- list(Parameters=MCMCRes$parameters,IRFs=IRFs,ModeParamTrans=dsgemode$par,ModeHessian=dsgemode$hessian,scalepar=scalepar,AcceptanceRate=MCMCRes$acceptRate,ObserveMat=ObserveMat,data=dsgedata,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds)
  class(dsgeret) <- "EDSGE"
  #
  return(dsgeret)
}

.DSGEParTransform <- function(parameters=NULL,parametersTrans=NULL,priorform,parbounds){
  if(class(parameters) == "numeric"){
    nparam <- as.numeric(length(parameters))
    parametersTrans <- numeric(nparam)
    #
    for(i in 1:nparam){
      if(priorform[i] == "Gamma" || priorform[i] == "IGamma"){
        if(class(parbounds[i,1]) != "numeric"){
          parametersTrans[i] <- log(parameters[i] - parbounds[i,1])
        }else{
          parametersTrans[i] <- log(parameters[i])
        }
      }else if(priorform[i] == "Beta"){
        parametersTrans[i] <- log((parameters[i] - parbounds[i,1])/(parbounds[i,2] - parameters[i]))
      }else{
        parametersTrans[i] <- parameters[i]
      }
    }
    #
    return(parametersTrans)
  }else{
    nparam <- as.numeric(length(parametersTrans))
    parameters <- numeric(nparam)
    #
    for(i in 1:nparam){
      if(priorform[i] == "Gamma" || priorform[i] == "IGamma"){
        if(class(parbounds[i,1]) != "numeric"){
          parameters[i] <- exp(parametersTrans[i]) + parbounds[i,1]
        }else{
          parameters[i] <- exp(parametersTrans[i])
        }
      }else if(priorform[i] == "Beta"){
        parameters[i] <- (parbounds[i,1] + parbounds[i,2]*exp(parametersTrans[i]))/(1 + exp(parametersTrans[i]))
      }else{
        parameters[i] <- parametersTrans[i]
      }
    }
    #
    return(parameters)
  }
}

.DSGEPriors <- function(parameters,parametersTrans,priorform,priorpars,parbounds,dsgelike){
  nparam <- as.numeric(length(parameters))
  dsgeposterior <- - dsgelike
  #
  for(i in 1:nparam){
    if(priorform[i] == "Normal"){
      dsgeposterior <- dsgeposterior + dnorm(parameters[i],mean = priorpars[i,1],sd = sqrt(priorpars[i,2]), log = TRUE)
    }else if(priorform[i] == "Gamma"){
      dsgeposterior <- dsgeposterior + dgamma(parameters[i],shape=priorpars[i,1],scale=priorpars[i,2],log=TRUE) + parametersTrans[i]
    }else if(priorform[i] == "IGamma"){
      dsgeposterior <- dsgeposterior + log( ( (priorpars[i,2]^priorpars[i,1])/gamma(priorpars[i,1]) )*( (parameters[i])^(-priorpars[i,1]-1) )*( exp(-priorpars[i,2]/parameters[i]) ) ) + parametersTrans[i]
    }else if(priorform[i] == "Beta"){
      if(is.na(parbounds[i,1]) || is.na(parbounds[i,2])){
        z = (parameters[i]-parbounds[i,1])/(parbounds[i,2]-parbounds[i,1])
        dsgeposterior <- dsgeposterior - log(parbounds[i,2] - parbounds[i,1]) - lbeta(priorpars[i,1],priorpars[i,2]) + ((priorpars[i,1]-1)*log(z)) + ((priorpars[i,2]-1)*log(1-z))
        dsgeposterior <- dsgeposterior + parametersTrans[i] - 2*log(1+exp(parametersTrans[i]))
      }else{
        z = (parameters[i]-parbounds[i,1])/(parbounds[i,2]-parbounds[i,1])
        dsgeposterior <- dsgeposterior - log(parbounds[i,2] - parbounds[i,1]) - lbeta(priorpars[i,1],priorpars[i,2]) + ((priorpars[i,1]-1)*log(z)) + ((priorpars[i,2]-1)*log(1-z))
        dsgeposterior <- dsgeposterior + log(parbounds[i,2] - parbounds[i,1]) + parametersTrans[i] - 2*log(1+exp(parametersTrans[i]))
      }
    }
  }
  #
  dsgeposterior <- - dsgeposterior
  #
  return(as.numeric(dsgeposterior))
}

.dsgeposteriorfn <- function(parametersTrans,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds){
  parameters <- .DSGEParTransform(NULL,parametersTrans,priorform,parbounds)
  dsgemats <- partomats(parameters)
  dsgesolved <- SDSGE(dsgemats$A,dsgemats$B,dsgemats$C,dsgemats$D,dsgemats$F,dsgemats$G,dsgemats$H,dsgemats$J,dsgemats$K,dsgemats$L,dsgemats$M,dsgemats$N)
  #
  StateMats <- .DSGEstatespace(dsgesolved$N,dsgesolved$P,dsgesolved$Q,dsgesolved$R,dsgesolved$S)
  #
  R <- matrix(0,nrow=ncol(ObserveMat),ncol=ncol(ObserveMat))
  dsgelike <- .Call("DSGEKalman", dsgedata,ObserveMat,StateMats$F,StateMats$G,dsgesolved$N,dsgemats$shocks,R,200, PACKAGE = "BMR", DUP = FALSE)
  dsgeposterior <- .DSGEPriors(parameters,parametersTrans,priorform,priorpars,parbounds,dsgelike$dsgelike)
  return(dsgeposterior)
}

.LaplaceMargLikelihood <- function(obj){
  PosteriorVal <- -obj$value
  MargLike <- PosteriorVal + (length(obj$par)*log(2*pi) + log(1/det(obj$hessian)))/2
  return(MargLike)
}

.DSGEMCMC <- function(dsgeopt,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds,parallel=FALSE){
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
    PrevLP <- (-1)*.dsgeposteriorfn(c(Draws[1,]),dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
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
    PropLP <- (-1)*.dsgeposteriorfn(c(proposal),dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
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

.DSGEMCMCMulti <- function(dsgeopt,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds,chains,cores){
  #
  cl <- makeCluster(cores)
  registerDoSNOW(cl)
  parallelsol <- 0
  parallelsol <- foreach(jj=1:chains, .packages=c("BMR")) %dopar% {
    .DSGEMCMC(dsgeopt,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds,TRUE)
  }
  #
  stopCluster(cl)
  #
  Draws <- matrix(NA,nrow=keep*chains,ncol=length(dsgeopt$par))
  accept <- numeric(chains)
  #
  for(j in 1:chains){
    Draws[((j-1)*keep+1):(j*keep),] <- parallelsol[[j]][[1]]
    accept[j] <- parallelsol[[j]][[2]]
  }
  #
  return=list(parameters=Draws,acceptRate=accept)
}

.MCMCDiagnostics <- function(parameters,chains){
  #
  keep <- nrow(parameters)/chains
  parArray <- array(0,dim=c(nrow(parameters),ncol(parameters),chains))
  for(i in 1:chains){
    parArray[,,i] <- parameters[((i-1)*keep+1):(i*keep),]
  }
  #
  meansAC <- numeric(dim(parArray)[2])
  for(i in 1:dim(parArray)[3]){
    for(j in 1:dim(parArray)[2]){
      meansAC[j] <- meansAC[j] + mean(parArray[,j,i])
    }
  }
  #
  meansAC <- meansAC/dim(parArray)[3]
  #
  B <- numeric(dim(parArray)[2])
  for(i in 1:dim(parArray)[3]){
    for(j in 1:dim(parArray)[2]){
      B[j] <- B[j] + (mean(parArray[,j,i]) - meansAC[j])^2
    }
  }
  #
  B <- dim(parArray)[1]/(dim(parArray)[3] - 1)*B
  #
  S <- numeric(dim(parArray)[2])
  for(i in 1:dim(parArray)[3]){
    for(j in 1:dim(parArray)[2]){
      S[j] <- S[j] + var(parArray[,j,i])
    }
  }
  W <- S/dim(parArray)[3]
  #
  #
  #
  rootR <- sqrt( ( ((dim(parArray)[1]-1)/dim(parArray)[1])*W + (1/dim(parArray)[1])*B ) / W)
  #
  return(rootR)
}