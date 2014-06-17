# 12/06/2014
BVARW.default <- function(mydata,cores,coefprior=NULL,p=4,constant=TRUE,irf.periods=20,keep=10000,burnin=1000,XiBeta=1,XiSigma=1,gamma=NULL){
  #
  kerr <- .bvarwerrors(mydata,cores,p,coefprior,constant,XiBeta,XiSigma,gamma)
  #
  kdata <- .bvarwdata(mydata,p,constant)
  #
  kprior <- .bvarwprior(kdata$Y,kdata$X,p,kdata$M,kdata$K,kdata$Yraw,kdata$Tp,kerr$coefprior,constant,kerr$XiBeta,kerr$XiSigma,kerr$gamma)
  #
  kreps <- 0
  if(cores==1){
    kreps <- .bvarwreplications(kdata$Y,kdata$X,kdata$Z,kprior$aPr,kprior$alpha,kprior$BVPr,kprior$vPr,kprior$SPr,kprior$Sigma,kdata$K,kdata$M,p,constant,kdata$Tp,keep,burnin,irf.periods)
  }else{
    kreps <- .bvarwreplicationsP(kdata$Y,kdata$X,kdata$Z,kprior$aPr,kprior$alpha,kprior$BVPr,kprior$vPr,kprior$SPr,kprior$Sigma,kdata$K,kdata$M,p,constant,kdata$Tp,keep,burnin,irf.periods,cores)
  }
  #
  bvarwret <- list(IRFs=kreps$IRFs,Beta=kreps$Beta,Sigma=kreps$Sigma,BDraws=kreps$BetaDraws,SDraws=kreps$SigmaDraws,data=mydata,constant=constant)
  class(bvarwret) <- "BVARW"
  return(bvarwret)
  #
}

.bvarwerrors <- function(mydata,cores,p,coefprior,constant,XiBeta,XiSigma,gamma){
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
  # Errors around user-given priors. First, Beta.
  #
  #No priors specified for the coefficients, set equal to random walk (in levels)
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
  # Prior on Beta Var
  if(class(XiBeta) == "numeric"){
    if(XiBeta <= 0){
      stop("XiBeta must be greater than zero.\n",call.=FALSE)
    }
  }else if(class(XiBeta) == "matrix"){
    if(constant==T){
      if(nrow(XiBeta) != ((ncol(mydata)*p+1)*ncol(mydata)) || ncol(XiBeta) != ((ncol(mydata)*p+1)*ncol(mydata))){
        stop("you have selected a full matrix for XiBeta.\n","Therefore, XiBeta must be of dimensions ",((ncol(mydata)*p+1)*ncol(mydata))," x ",((ncol(mydata)*p+1)*ncol(mydata)),".\n",call.=FALSE)
      }
    }else{
      if(nrow(XiBeta) != ((ncol(mydata)*p)*ncol(mydata)) || ncol(XiBeta) != ((ncol(mydata)*p)*ncol(mydata))){
        stop("you have selected a full matrix for XiBeta.\n","Therefore, XiBeta must be of dimensions ",((ncol(mydata)*p)*ncol(mydata))," x ",((ncol(mydata)*p)*ncol(mydata)),".\n",call.=FALSE)
      }
    }
  }else{
    if(constant==T){
      stop("unrecognised form for XiBeta.\n","XiBeta must be a ",((ncol(mydata)*p+1)*ncol(mydata))," x ",((ncol(mydata)*p+1)*ncol(mydata))," matrix, or single numeric value.\n",call.=FALSE)
    }else{
      stop("unrecognised form for XiBeta.\n","XiBeta must be a ",((ncol(mydata)*p)*ncol(mydata))," x ",((ncol(mydata)*p)*ncol(mydata))," matrix, or single numeric value.\n",call.=FALSE)
    }
  }
  #
  # Now Sigma
  # 
  if(class(XiSigma) == "numeric"){
    if(XiSigma <= 0){
      stop("XiSigma must be greater than zero.\n",call.=FALSE)
    }
  }else if(class(XiSigma) == "matrix"){
    if(nrow(XiSigma) != (ncol(mydata)) || ncol(XiSigma) != (ncol(mydata))){
      stop("you have selected a full matrix for XiSigma.\n","Therefore, XiSigma must be of dimensions ",(ncol(mydata))," x ",(ncol(mydata)),".\n",call.=FALSE)
    }
  }else{
    stop("unrecognised form for XiSigma.\n","XiSigma must be a ",(ncol(mydata))," x ",(ncol(mydata))," matrix, or single numeric value.\n",call.=FALSE)
  }
  #
  # Prior COV Degrees of Freedom
  if(class(gamma)=="NULL"){
    gamma <- ncol(mydata) + 1
  }else if(gamma <= ncol(mydata)){
    stop("with this data, the minimum value for gamma is ",(ncol(mydata)+1),".\n",call.=FALSE)
  }else{gamma <- gamma}
  #
  #
  #
  return=list(coefprior=coefprior,XiBeta=XiBeta,XiSigma=XiSigma,gamma=gamma)
}

.bvarwdata <- function(mydata,p,constant){
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
  Z<-kronecker(diag(M),X)
  Y <- Yraw[(p+1):nrow(Yraw),]
  #
  return=list(Y=Y,X=X,Z=Z,M=M,K=K,Yraw=Yraw,Tp=Tp)
}

.bvarwprior <- function(Y,X,p,M,K,Yraw,Tp,coefprior,constant,XiBeta,XiSigma,gamma){
  Beta <- solve.qr(qr(as.matrix(X)),as.matrix(Y))
  alpha <- c(Beta)
  #
  Epsilon <- as.matrix(Y) - as.matrix(X)%*%Beta
  Sigma <- (t(Epsilon)%*%Epsilon)/(Tp - K)
  #
  #
  #
  BPr <- 0
  if(class(coefprior)=="numeric"){
    if(constant==T){
      BPr<-rbind(rep(0,M),diag(M),matrix(rep(0,(p-1)*M*M),nrow=(p-1)*M,ncol=M))
    }else{
      BPr<-rbind(diag(M),matrix(rep(0,(p-1)*M*M),nrow=(p-1)*M,ncol=M))
    }
    for(i in 1:M){BPr[(i+1),i]<-coefprior[i]}
  }else{BPr <- coefprior}
  aPr <- c(BPr)
  #
  BVPr<-0
  if(class(XiBeta)=="numeric"){
    BVPr <- diag(K*M)*XiBeta
  }else{BVPr <- XiBeta}
  #
  #
  #
  SPr<-0
  if(class(XiSigma)=="numeric"){
    SPr <- diag(M)*XiSigma
  }else{SPr <- XiSigma}
  #
  vPr <- gamma
  #
  return=list(aPr=aPr,alpha=alpha,BVPr=BVPr,vPr=vPr,SPr=SPr,Sigma=Sigma)
}

.bvarwreplications <- function(Y,X,Z,aPr,alpha,BVPr,vPr,SPr,Sigma,K,M,p,constant,Tp,keep,burnin,irf.periods){
  #
  ImpStore <- 0
  #
  cat('Starting Gibbs C++, ', date(),'. \n', sep="")
  #RepsRun <- WBVARReps(Sigma,as.matrix(X),as.matrix(Z),as.matrix(Y),matrix(aPr,ncol=1),SPr,vPr,BVPr,Tp,M,K,burnin,keep)
  RepsRun <- .Call("WBVARReps", Sigma,as.matrix(X),as.matrix(Z),as.matrix(Y),matrix(aPr,ncol=1),SPr,vPr,BVPr,Tp,M,K,burnin,keep, PACKAGE = "BMR", DUP = FALSE)
  cat('C++ reps finished, ', date(),'. Now generating IRFs. \n', sep="")
  #
  kcons <- 0; if(constant==T){kcons<-1}
  #ImpStore <- WBVARIRFs(M,K,kcons,keep,irf.periods,RepsRun$Beta,RepsRun$Sigma)
  ImpStore <- .Call("WBVARIRFs", M,K,kcons,keep,irf.periods,RepsRun$Beta,RepsRun$Sigma, PACKAGE = "BMR", DUP = FALSE)
  ImpStore <- ImpStore$ImpStore
  IRFStore <- array(NA,dim=c(M,M,irf.periods,keep))
  for(i in 1:keep){
    IRFStore[,,,i] <- ImpStore[,,((i-1)*irf.periods+1):(i*irf.periods)]
  }
  #
  BetaMean <- apply(RepsRun$Beta,c(1,2),mean)
  SigmaMean <- apply(RepsRun$Sigma,c(1,2),mean)
  #
  #
  IRFSorted<-apply(IRFStore,c(3,1,2),sort)
  #
  IRFSorted<-aperm(IRFSorted,c(2,3,1,4))
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
  return=list(Beta=BetaMean,Sigma=SigmaMean,BetaDraws=RepsRun$Beta,SigmaDraws=RepsRun$Sigma,IRFs=IRFSorted)
}

.bvarwreplicationsP <- function(Y,X,Z,aPr,alpha,BVPr,vPr,SPr,Sigma,K,M,p,constant,Tp,keep,burnin,irf.periods,NCore){
  #
  require(doSNOW)
  #
  ImpStore <- 0
  #
  cat('Starting Gibbs C++, ', date(),'. \n', sep="")
  #
  #RepsRunB <- WBVARRepsB(Sigma,as.matrix(X),as.matrix(Z),as.matrix(Y),matrix(aPr,ncol=1),SPr,vPr,BVPr,Tp,M,K,burnin)
  RepsRunB <- .Call("WBVARRepsB", Sigma,as.matrix(X),as.matrix(Z),as.matrix(Y),matrix(aPr,ncol=1),SPr,vPr,BVPr,Tp,M,K,burnin, PACKAGE = "BMR", DUP = FALSE)
  #
  cl <- makeCluster(NCore)
  registerDoSNOW(cl)
  keeppar <- ceiling(keep/NCore)
  solpar <- 0
  solpar <- foreach(jj=1:NCore, .packages=c("BMR")) %dopar% {
    .RepsBFn(RepsRunB$Sigma,as.matrix(X),as.matrix(Z),as.matrix(Y),matrix(aPr,ncol=1),SPr,vPr,BVPr,Tp,M,K,keeppar)
  }
  #
  stopCluster(cl)
  #
  BetaArray <- array(0,dim=c(K,M,keeppar*NCore))
  SigmaArray <- array(0,dim=c(M,M,keeppar*NCore))
  for(j in 1:NCore){
    BetaArray[,,((j-1)*keeppar+1):(j*keeppar)] <- solpar[[j]][[1]]
    SigmaArray[,,((j-1)*keeppar+1):(j*keeppar)] <- solpar[[j]][[2]]
  }
  BetaArray <- BetaArray[,,1:keep]
  SigmaArray <- SigmaArray[,,1:keep]
  #
  cat('C++ reps finished, ', date(),'. Now generating IRFs. \n', sep="")
  #
  kcons <- 0; if(constant==T){kcons<-1}
  ImpStore <- .Call("WBVARIRFs", M,K,kcons,keep,irf.periods,BetaArray,SigmaArray, PACKAGE = "BMR", DUP = FALSE)
  ImpStore <- ImpStore$ImpStore
  IRFStore <- array(NA,dim=c(M,M,irf.periods,keep))
  for(i in 1:keep){
    IRFStore[,,,i] <- ImpStore[,,((i-1)*irf.periods+1):(i*irf.periods)]
  }
  #
  BetaMean <- apply(BetaArray,c(1,2),mean)
  SigmaMean <- apply(SigmaArray,c(1,2),mean)
  #
  #
  IRFSorted<-apply(IRFStore,c(3,1,2),sort)
  #
  IRFSorted<-aperm(IRFSorted,c(2,3,1,4))
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
  return=list(Beta=BetaMean,Sigma=SigmaMean,BetaDraws=BetaArray,SigmaDraws=SigmaArray,IRFs=IRFSorted)
}

.RepsBFn <- function(Sigma,X,Z,Y,aPr,SPr,vPr,BVPr,Tp,M,K,keeppar){
  #
  #Res <- WBVARRepsK(Sigma,X,Z,Y,aPr,SPr,vPr,BVPr,Tp,M,K,keeppar)
  Res <- .Call("WBVARRepsK", Sigma,X,Z,Y,aPr,SPr,vPr,BVPr,Tp,M,K,keeppar, PACKAGE = "BMR", DUP = FALSE)
  #
  return(list(Res$Beta,Res$Sigma))
}