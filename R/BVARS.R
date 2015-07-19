################################################################################
##
##   R package BMR by Keith O'Hara Copyright (C) 2011, 2012, 2013, 2014, 2015
##   This file is part of the R package BMR.
##
##   The R package BMR is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 2 of the License, or
##   (at your option) any later version.
##
##   The R package BMR is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
################################################################################

# 07/20/2015

BVARS.default <- function(mydata,psiprior=NULL,coefprior=NULL,p=4,irf.periods=20,keep=10000,burnin=1000,XiPsi=1,HP1=0.5,HP4=2,gamma=NULL){
  #
  kerr <- .bvarserrors(mydata,p,coefprior,psiprior,XiPsi,HP1,HP4,gamma)
  #
  kdata <- .bvarsdata(mydata,p)
  #
  kprior <- .bvarsprior(kdata$Y,kdata$X,kdata$d,kdata$dX,p,kdata$M,kdata$K,kdata$Yr,kdata$Tr,kdata$Tp,kerr$coefprior,kerr$psiprior,XiPsi,HP1,HP4)
  #
  kreps <- .bvarsreplications(kdata$Y,kdata$X,kdata$d,kdata$dX,kprior$yd,kprior$Zd,kprior$PsiPr,kprior$invPsiVPr,kprior$BPr,kprior$Beta,kprior$invBVPr,kprior$Sigma,kprior$SigmaML,kerr$gamma,kdata$K,kdata$M,p,kdata$Tp,keep,burnin,irf.periods)
  #
  bvarsret <- list(IRFs=kreps$IRFs,Psi=kreps$Psi,Beta=kreps$Beta,Sigma=kreps$Sigma,PDraws=kreps$PsiDraws,BDraws=kreps$BetaDraws,SDraws=kreps$SigmaDraws,data=mydata)
  class(bvarsret) <- "BVARS"
  return(bvarsret)
  #
}

.bvarserrors <- function(mydata,p,coefprior,psiprior,XiPsi,HP1,HP4,gamma){
  if(ncol(mydata)<2){
    stop("need more than 1 variable in your data.\n",call.=FALSE)
  }
  if(p<1){
    stop("need at least 1 lag; set p higher.\n",call.=FALSE)
  }
  if(p>nrow(mydata)){
    stop("need more data points than lags.\n",call.=FALSE)
  }
  for(i in 1:ncol(mydata)){
    if(sum(is.na(mydata[,i]))){
      stop("no missing observations allowed.\n",call.=FALSE)
    }
  }
  #
  # Errors around user-given priors. Psi first.
  #
  if(class(psiprior)=="NULL"){
    stop("a prior on Psi MUST be specified.\n",call.=FALSE)
  }
  #
  # Errors around user-given priors. Psi first.
  #
  if(class(psiprior)=="matrix"){
    if(length(c(psiprior)) != (ncol(mydata))){
      stop("the matrix of priors on psi is not of the appropriate dimensions, which are ",1," x ",ncol(mydata),".\n",call.=FALSE)
    }
    psiprior <- c(psiprior)
  }else if(class(psiprior)=="numeric"){
    if(length(psiprior) != ncol(mydata)){
      stop("the prior on psi is not of the appropriate dimensions; please provide a numeric vector of length ",ncol(mydata),".\n",call.=FALSE)
    }
  }else{psiprior <- psiprior}
  #
  # 
  if(class(XiPsi) == "numeric"){
    if(XiPsi <= 0){
      stop("XiPsi must be greater than zero.\n",call.=FALSE)
    }
  }else if(class(XiPsi) == "matrix"){
    if(nrow(XiPsi) != (ncol(mydata)) || ncol(XiPsi) != (ncol(mydata))){
      stop("you have selected a full matrix for XiPsi.\n","Therefore, XiPsi must be of dimensions ",(ncol(mydata))," x ",(ncol(mydata)),".\n",call.=FALSE)
    }
  }else{
    stop("unrecognised form for XiPsi.\n","XiPsi must be a ",(ncol(mydata))," x ",(ncol(mydata))," matrix, or single numeric value.\n",call.=FALSE)
  }
  #
  # Now Beta
  #No priors specified for the coefficients, set equal to random walk (in levels)
  if(class(coefprior)=="NULL"){
    coefprior <- c(rep(1,ncol(mydata)))
  }
  if(class(coefprior)=="matrix"){
    if(length(c(coefprior)) != ((ncol(mydata)*p)*ncol(mydata))){
      stop("you have opted to give a full prior matrix on Beta. However,\n", "it is not of the appropriate dimensions, which are ",(ncol(mydata)*p)," x ",ncol(mydata),".\n",call.=FALSE)
    }
  }
  if(class(coefprior)=="numeric"){
    if(length(coefprior) != ncol(mydata)){
      stop("you have opted to give a prior on the first-own lags in Beta. However,\n", "it is not of the appropriate dimensions; please provide a numeric vector of length ",ncol(mydata),".\n",call.=FALSE)
    }
  }else{coefprior <- coefprior}
  #
  if(HP4 <= 0){
    stop("HP4 must be greater than zero.\n",call.=FALSE)
  }
  #
  # Sigma
  # Prior COV Degrees of Freedom
  if(class(gamma)=="NULL"){
    gamma <- ncol(mydata) + 1
  }else if(gamma <= ncol(mydata)){
    stop("with this data, the minimum value for gamma is ",(ncol(mydata)+1),".\n",call.=FALSE)
  }else{gamma <- gamma}
  #
  #
  #
  return=list(coefprior=coefprior,psiprior=psiprior,gamma=gamma)
}

.bvarsdata <- function(mdata,p){
  Tr <- as.numeric(dim(mdata)[1])
  M <- as.numeric(dim(mdata)[2])
  #
  Yr <- as.matrix(mdata,ncol=M)
  #
  X <- embed(Yr,p+1); X <- X[,(M+1):ncol(X)]
  #
  K <- as.numeric(dim(X)[2])
  Y <- Yr[(p+1):nrow(Yr),]
  Tp <- Tr - p
  #
  d <- matrix(rep(1,Tp))
  dX <- matrix(rep(-1,Tp*p),ncol=p)
  #
  return=list(Y=Y,X=X,d=d,dX=dX,M=M,K=K,Yr=Yr,Tr=Tr,Tp=Tp)
}

.bvarsprior <- function(Y,X,d,dX,p,M,K,Yr,Tr,Tp,coefprior,psiprior,XiPsi,HP1,HP4){
  #
  # Psi
  #
  PsiPr <- matrix(psiprior,ncol=M)
  #
  if(class(XiPsi)=="numeric"){
    PsiVPr <- diag(M)*XiPsi
  }else{PsiVPr <- XiPsi}
  invPsiVPr <- solve(PsiVPr)
  #
  Psi.int <- colMeans(Y); Psi.int <- matrix(Psi.int,ncol=M)
  #
  # Beta
  #
  if(class(coefprior)=="numeric"){
    BPr<-rbind(diag(M),matrix(rep(0,(p-1)*M*M),nrow=(p-1)*M,ncol=M))
    for(i in 1:M){BPr[i,i]<-coefprior[i]}
  }else{BPr <- coefprior}
  #
  BVPr <- matrix(0,M*p,M*p)
  for(i in 1:p){
    BVPr[(1 + (M*(i-1))):(M*i),(1 + (M*(i-1))):(M*i)] <- (HP1/(i^(HP4)))*diag(M)
  }
  invBVPr <- solve(BVPr)
  #
  # Sigma
  #
  yd <- Y - d%*%Psi.int
  Zd <- X + dX%*%kronecker(diag(p),Psi.int)
  Beta <- solve.qr(qr(Zd),yd)
  #
  Epsilon <- yd - Zd%*%Beta
  Sigma <- (t(Epsilon)%*%Epsilon)/(Tp + M + 1)
  SigmaML <- Sigma
  #
  return=list(yd=yd,Zd=Zd,PsiPr=PsiPr,invPsiVPr=invPsiVPr,BPr=BPr,Beta=Beta,invBVPr=invBVPr,Sigma=Sigma,SigmaML=SigmaML)
}

.bvarsreplications <- function(Y,X,d,dX,yd,Zd,PsiPr,invPsiVPr,BPr,Beta,invBVPr,Sigma,SigmaML,gamma,K,M,p,Tp,keep,burnin,irf.periods){
  #
  ImpStore <- 0
  #
  message('Starting Gibbs C++, ', date(),'.', sep="")
  RepsRun <- .Call("SBVARReps", as.matrix(X),as.matrix(Y),d,dX,yd,Zd,PsiPr,invPsiVPr,BPr,Beta,invBVPr,Sigma,SigmaML,gamma,Tp,M,p,burnin,keep, PACKAGE = "BMR")
  message('C++ reps finished, ', date(),'. Now generating IRFs.', sep="")
  #
  ImpStore <- .Call("SBVARIRFs", M,K,keep,irf.periods,RepsRun$Beta,RepsRun$Sigma, PACKAGE = "BMR")
  #
  ImpStore <- ImpStore$ImpStore
  ImpStore2 <- array(NA,dim=c(M,M,irf.periods,keep))
  for(i in 1:keep){
    ImpStore2[,,,i] <- ImpStore[,,((i-1)*irf.periods+1):(i*irf.periods)]
  }
  #
  PsiMean <- apply(RepsRun$Psi,c(1,2),mean)
  BetaMean <- apply(RepsRun$Beta,c(1,2),mean)
  SigmaMean <- apply(RepsRun$Sigma,c(1,2),mean)
  #
  #
  ImpSorted <- apply(ImpStore2,c(3,1,2),sort)
  #
  ImpSorted <- aperm(ImpSorted,c(2,3,1,4))
  #
  rownames(PsiMean) <- "Psi"
  #
  BNames <- character(length=nrow(BetaMean))
  kcon <- 1
  for(i in 1:p){  
    for(j in 1:M){
      BNames[kcon] <- paste("Eq ",j,", lag ",i,sep="")
      kcon <- kcon + 1
    }
  }
  rownames(BetaMean) <- BNames
  if(class(colnames(Y)) == "character"){
    colnames(BetaMean) <- colnames(Y)
    colnames(PsiMean) <- colnames(Y)
  }
  #
  return=list(Psi=PsiMean,Beta=BetaMean,Sigma=SigmaMean,PsiDraws=RepsRun$Psi,BetaDraws=RepsRun$Beta,SigmaDraws=RepsRun$Sigma,IRFs=ImpSorted)
}