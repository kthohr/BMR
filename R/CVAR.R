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

CVAR.default <- function(mydata,p=4,constant=TRUE,irf.periods=20,boot=10000){
  #
  kerr <- .cvarerrors(mydata,p)
  #
  kdata <- .cvardata(mydata,p,constant)
  #
  kprior <- .cvarsetup(kdata$Y,kdata$X,kdata$K,kdata$Tp)
  #
  kreps <- .cvarboot(kdata$Y,kdata$X,kprior$Beta,kprior$Sigma,kprior$Resids,kdata$K,kdata$M,p,constant,kdata$Tp,boot,irf.periods)
  #
  cvarret <- list(IRFs=kreps$IRFs,Beta=kreps$Beta,Sigma=kreps$Sigma,BDraws=kreps$BetaDraws,SDraws=kreps$SigmaDraws,data=mydata,constant=constant)
  class(cvarret) <- "CVAR"
  return(cvarret)
  #
}

.cvarerrors <- function(mydata,p){
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
}

.cvardata <- function(mydata,p,constant){
  Tr <- as.numeric(dim(mydata)[1])
  M <- as.numeric(dim(mydata)[2])
  #
  Yr <- as.matrix(mydata,ncol=M)
  #
  X <- embed(Yr,p+1); X <- X[,(M+1):ncol(X)]
  if(constant == TRUE){X<-cbind(rep(1,(Tr-p)),X)}
  #
  K <- as.numeric(dim(X)[2])
  Y <- Yr[(p+1):nrow(Yr),]
  Tp <- Tr - p
  #
  return=list(Y=Y,X=X,M=M,K=K,Tp=Tp)
}

.cvarsetup <- function(Y,X,K,Tp){
  Beta <- solve.qr(qr(as.matrix(X)),as.matrix(Y))
  #
  Resids <- as.matrix(Y) - as.matrix(X)%*%Beta
  #
  Sigma<-(1/(Tp-K))*(t(Resids)%*%Resids)
  #
  return=list(Beta=Beta,Sigma=Sigma,Resids=Resids)
}

.cvarboot <- function(Y,X,Beta,Sigma,Resids,K,M,p,constant,Tp,boot,irf.periods){
  #
  ImpStore <- 0
  ResidDraws <- array(0,dim=c(Tp,M,boot))
  #
  for(i in 1:boot){
    SPoints <- sample(c(1:(Tp)),replace=TRUE)
    ResidDraws[,,i] <- Resids[SPoints,]
  }
  #
  kcons <- 0; if(constant==T){kcons<-1}
  #
  message('Starting C++, ', date(),'.', sep="")
  RepsRun <- .Call("CVARReps", Beta,Sigma,as.matrix(X),as.matrix(Y),Tp,M,K,kcons,boot,ResidDraws, PACKAGE = "BMR")
  message('C++ reps finished, ', date(),'. Now getting IRFs.', sep="")
  #
  ImpStore <- .Call("CVARIRFs", M,K,kcons,boot,irf.periods,RepsRun$Beta,RepsRun$Sigma, PACKAGE = "BMR")
  ImpStore <- ImpStore$ImpStore
  ImpStore2 <- array(NA,dim=c(M,M,irf.periods,boot))
  for(i in 1:boot){
    ImpStore2[,,,i] <- ImpStore[,,((i-1)*irf.periods+1):(i*irf.periods)]
  }
  #
  BNames <- character(length=nrow(Beta))
  kcon <- 1; if(constant==T){BNames[1] <- "Constant"};if(constant==T){kcon<-2}
  for(i in 1:p){  
    for(j in 1:M){
      BNames[kcon] <- paste("Eq ",j,", lag ",i,sep="")
      kcon <- kcon + 1
    }
  }
  rownames(Beta) <- BNames
  if(class(colnames(Y)) == "character"){
    colnames(Beta) <- colnames(Y)
  }
  #
  ImpSorted <- apply(ImpStore2,c(3,1,2),sort)
  #
  ImpSorted <- aperm(ImpSorted,c(2,3,1,4))
  #
  return=list(Beta=Beta,Sigma=Sigma,BetaDraws=RepsRun$Beta,SigmaDraws=RepsRun$Sigma,IRFs=ImpSorted)
}