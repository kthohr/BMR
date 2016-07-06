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

BVARM.default <- function(mydata,coefprior=NULL,p=4,constant=TRUE,irf.periods=20,keep=10000,burnin=1000,VType=1,decay="H",HP1=0.5,HP2=0.5,HP3=1,HP4=2){
  #
  kerr <- .bvarmerrors(mydata,p,coefprior,constant,VType,decay,HP4)
  #
  kdata <- .bvarmdata(mydata,p,constant)
  #
  kprior <- .bvarmprior(kdata$Y,kdata$X,p,kdata$M,kdata$K,kdata$Yraw,kdata$Tp,kerr$coefprior,constant,VType,decay,HP1,HP2,HP3,HP4)
  #
  kreps <- .bvarmreplications(kdata$Y,kdata$Z,kprior$aPr,kprior$alpha,kprior$BVPr,kprior$Sigma,kdata$K,kdata$M,p,constant,keep,burnin,irf.periods)
  #
  bvarmret <- list(IRFs=kreps$IRFs,BetaVPr=kprior$BVPr,Beta=kreps$Beta,Sigma=kreps$Sigma,BDraws=kreps$BetaDraws,data=mydata,constant=constant)
  class(bvarmret) <- "BVARM"
  return(bvarmret)
  #
}

.bvarmerrors <- function(mydata,p,coefprior,constant,VType,decay,HP4){
  #
  if(ncol(mydata)<2){
    stop("need more than 1 variable.\n",call.=FALSE)
  }
  if(p<1){
    stop("need at least 1 lag.\n",call.=FALSE)
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
  # Errors around prior coefs
  #
  # If no prior is provided for the coefficients, set to a random walk (in levels)
  if(class(coefprior)=="NULL"){
    coefprior <- c(rep(1,ncol(mydata)))
  }
  if(class(coefprior)=="matrix"){
    if(constant==T){
      if(length(c(coefprior)) != ((ncol(mydata)*p + 1)*ncol(mydata))){
        stop("you have opted to give a full prior matrix on Beta. However,\n", "it is not of the appropriate dimensions, which are ",(ncol(mydata)*p + 1)," x ",ncol(mydata),".\n",call.=FALSE)
      }
    }else{
      if(length(c(coefprior)) != ((ncol(mydata)*p)*ncol(mydata))){
        stop("you have opted to give a full prior matrix on Beta. However,\n", "it is not of the appropriate dimensions, which are ",(ncol(mydata)*p)," x ",ncol(mydata),".\n",call.=FALSE)
      }
    }
  }
  if(class(coefprior)=="numeric"){
    if(length(coefprior) != ncol(mydata)){
      stop("you have opted to give a prior on the first-own lags in Beta. However,\n", "it is not of the appropriate dimensions; please provide a numeric vector of length ",ncol(mydata),".\n",call.=FALSE)
    }
  }else{coefprior <- coefprior}
  #
  #
  # Some possible errors around hyperparameters and such
  if(VType != 1 && VType != 2){
    stop("VType can only be 1 or 2.\n",call.=FALSE)
  }
  if(decay != "H" && decay != "G"){
    stop("decay must be H or G (and quotation marks both sides!).\n",call.=FALSE)
  }
  if(HP4 <= 0){
    stop("HP4 must be greater than zero.\n",call.=FALSE)
  }
  #
  return=list(coefprior=coefprior)
}

.bvarmdata <- function(mydata,p,constant){
  Tr <- as.numeric(dim(mydata)[1])
  M <- as.numeric(dim(mydata)[2])
  #
  Yraw <- as.matrix(mydata,ncol=M)
  #
  X <- embed(Yraw,p+1); X <- X[,(M+1):ncol(X)]
  if(constant == TRUE){X<-cbind(rep(1,(Tr-p)),X)}
  #
  K <- as.numeric(dim(X)[2])
  Z <- kronecker(diag(M),X)
  Y <- Yraw[(p+1):nrow(Yraw),]
  Tp <- Tr - p
  #
  return=list(Y=Y,X=X,Z=Z,M=M,K=K,Yraw=Yraw,Tp=Tp)
}

.bvarmprior <- function(Y,X,p,M,K,Yraw,Tp,coefprior,constant,VType,decay,HP1,HP2,HP3,HP4){
  Beta <- solve.qr(qr(as.matrix(X)),as.matrix(Y))
  alpha <- c(Beta)
  #
  if(class(coefprior)=="numeric"){
    if(constant==T){
      BPr<-rbind(rep(0,M),diag(M),matrix(rep(0,(p-1)*M*M),nrow=(p-1)*M,ncol=M))
      for(i in 1:M){BPr[(i+1),i]<-coefprior[i]}
    }else{
      BPr<-rbind(diag(M),matrix(rep(0,(p-1)*M*M),nrow=(p-1)*M,ncol=M))
      for(i in 1:M){BPr[i,i]<-coefprior[i]}
    }
  }else{BPr <- coefprior}
  aPr <- c(BPr)
  #
  Sigma<-rep(0,M)
  #
  for(i in 1:M){
    XAR <- embed(Yraw[,i],p+1)
    YAR <- matrix(XAR[,1])
    XAR <- XAR[,2:ncol(XAR)]; XAR <- as.matrix(XAR,ncol=p); 
    if(constant==T){ XAR <- cbind(as.matrix(rep(1,nrow(YAR)),ncol=1),XAR) }
    #
    alphaAR <- solve.qr(qr(XAR),YAR)
    if(constant==T){
      Sigma[i] <- (1/(Tp-p-1))*t(YAR - XAR%*%alphaAR)%*%(YAR - XAR%*%alphaAR)
    }else{
      Sigma[i] <- (1/(Tp-p))*t(YAR - XAR%*%alphaAR)%*%(YAR - XAR%*%alphaAR)
    }
  }
  #
  # Prior covariance matrix of alpha
  #
  BVPr <- matrix(0,nrow=K,ncol=M)
  #
  if(VType==1){
    Cons <- 0
    if(constant==TRUE){
      BVPr[1,] <- Sigma*HP3
      Cons <- 1
    }
    #
    for(i in 1:p){
      #j is the column, k is the row
      for(j in 1:M){
        for(k in 1:M){
          BVPr[((i-1)*M + k + Cons),j] <- (HP2*Sigma[j])/((i^2)*Sigma[k])
        }
        BVPr[((i-1)*M + j + Cons),j] <- HP1/(i^2)
      }
    }
  }else{
    Cons <- 0
    if(constant==TRUE){
      BVPr[1,] <- HP1*HP3
      Cons <- 1
    }
    # Lag decay function
    if(decay=="H"){
      LDecay <- function(x,HP4){
        VDecay <- x^(HP4)
        return(VDecay)
      }
    }else{
      LDecay <- function(x,HP4){
        VDecay <- HP4^((-x) + 1)
        return(VDecay)
      }
    }
    #
    for(i in 1:p){
      #j is the column, k is the row
      for(j in 1:M){
        for(k in 1:M){
          BVPr[((i-1)*M + k + Cons),j] <- HP1*(HP2*Sigma[k])/((LDecay(i,HP4))*Sigma[j])
        }
        BVPr[((i-1)*M + j + Cons),j] <- HP1/(LDecay(i,HP4))
      }
    }
  }
  #
  #
  #
  BVPr <- c(BVPr)
  BVPr <- diag(BVPr)
  #
  Sigma <- c(Sigma)
  Sigma <- diag(Sigma)
  #
  return=list(aPr=aPr,alpha=alpha,BVPr=BVPr,Sigma=Sigma)
}

.bvarmreplications <- function(Y,Z,aPr,alpha,BVPr,Sigma,K,M,p,constant,keep,burnin,irf.periods){
  #
  ImpStore <- 0
  #
  SigmaO<-Sigma
  shock<-chol(Sigma)
  #
  Sigma <- kronecker(solve(Sigma),diag(nrow(Y)))
  #
  message('Starting Gibbs C++, ', date(),'.', sep="")
  RepsRun <- .Call("MBVARReps", Sigma,as.matrix(Z),as.matrix(Y),matrix(aPr,ncol=1),BVPr,M,K,burnin,keep, PACKAGE = "BMR")
  message('C++ reps finished, ', date(),'. Now generating IRFs.', sep="")
  #
  kcons <- 0; if(constant==T){kcons<-1}
  ImpStore <- .Call("MBVARIRFs", shock,M,K,kcons,keep,irf.periods,RepsRun$Beta, PACKAGE = "BMR")
  ImpStore <- ImpStore$ImpStore
  ImpStore2 <- array(NA,dim=c(M,M,irf.periods,keep))
  for(i in 1:keep){
    ImpStore2[,,,i] <- ImpStore[,,((i-1)*irf.periods+1):(i*irf.periods)]
  }
  #
  BetaMean <- apply(RepsRun$Beta,c(1,2),mean)
  #
  BNames <- character(length=nrow(BetaMean))
  kcon <- 1; if(constant==T){BNames[1] <- "Constant"};if(constant==T){kcon<-2}
  for(i in 1:p){  
    for(j in 1:M){
      BNames[kcon] <- paste("Eq ",j,", lag ",i,sep="")
      kcon <- kcon + 1
    }
  }
  rownames(BetaMean) <- BNames
  if(class(colnames(Y)) == "character"){
    colnames(BetaMean) <- colnames(Y)
  }
  #
  ImpSorted <- apply(ImpStore2,c(3,1,2),sort)
  #
  ImpSorted <- aperm(ImpSorted,c(2,3,1,4))
  #
  return=list(Beta=BetaMean,BetaDraws=RepsRun$Beta,Sigma=SigmaO,IRFs=ImpSorted)
}