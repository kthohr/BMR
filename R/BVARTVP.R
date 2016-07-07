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

BVARTVP.default <- function(mydata,timelab=NULL,coefprior=NULL,tau=NULL,p=4,irf.periods=20,irf.points=NULL,keep=10000,burnin=5000,XiBeta=1,XiQ=0.01,gammaQ=NULL,XiSigma=1,gammaS=NULL){
    #
    kerr <- .bvartvperror(mydata,p,coefprior,tau,XiBeta,XiQ,gammaQ,XiSigma,gammaS)
    #
    kdata <- .bvartvpdata(mydata,p,tau,timelab,irf.points)
    #
    kprior <- .bvartvpprior(kdata$Y,kdata$Z,p,kdata$M,kdata$K,kdata$kT,kerr$coefprior,tau,keep,irf.periods,kerr$XiBeta,kerr$XiQ,kerr$gammaQ,kerr$XiSigma,kerr$gammaS)
    #
    message('Finished Prior.')
    #
    kreps <- .bvartvpreps(kdata$Y,kdata$y,kdata$Z,kprior$B0Pr,kprior$B0VPr,kprior$invB0VPr,
                          kprior$QPr,kprior$QVPr,kprior$SPr,kprior$SVPr,kprior$QDraw,kprior$Qchol,
                          kprior$SDraw,kprior$invSDraw,kdata$K,kdata$M,p,kdata$kT,
                          keep,burnin,kdata$timelab,kdata$nIRFs,kdata$irf.points,irf.periods)
    #
    bvartvpret <- list(IRFs=kreps$IRFs,Beta=kreps$BetaMean,Q=kreps$QMean,Sigma=kreps$SigmaMean,BDraws=kreps$Betas,QDraws=kreps$QDraws,SDraws=kreps$SDraws,data=mydata,irf.points=irf.points,tau=tau)
    class(bvartvpret) <- "BVARTVP"
    return(bvartvpret)
}

.bvartvperror <- function(mydata,p,coefprior,tau,XiBeta,XiQ,gammaQ,XiSigma,gammaS){
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
    # Errors around user-given priors. First, Beta.
    #
    # Prior on Beta
    if(class(coefprior)=="NULL"){
        coefprior <- c(rep(1,ncol(mydata)))
    }
    if(class(coefprior)=="matrix"){
        if(length(c(coefprior)) != ((ncol(mydata)*p + 1)*ncol(mydata))){
            stop("you have opted to give a full prior matrix on Beta. However,\n", "it is not of the appropriate dimensions, which are ",(ncol(mydata)*p + 1)," x ",ncol(mydata),".\n",call.=FALSE)
        }
    }
    if(class(coefprior)=="numeric"){
        if(length(coefprior) != ncol(mydata) && class(tau) == "NULL"){
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
        if(nrow(XiBeta) != ((ncol(mydata)*p+1)*ncol(mydata)) || ncol(XiBeta) != ((ncol(mydata)*p+1)*ncol(mydata))){
            stop("you have selected a full matrix for XiBeta.\n","Therefore, XiBeta must be of dimensions ",((ncol(mydata)*p+1)*ncol(mydata))," x ",((ncol(mydata)*p+1)*ncol(mydata)),".\n",call.=FALSE)
        }
    }else{
        stop("unrecognised form for XiBeta.\n","XiBeta must be a ",((ncol(mydata)*p+1)*ncol(mydata))," x ",((ncol(mydata)*p+1)*ncol(mydata))," matrix, or single numeric value.\n",call.=FALSE)
    }
    #
    # Now Q
    # 
    if(class(XiQ) == "numeric"){
        if(XiQ > 0){
            XiQ <- XiQ
        }else if(XiQ <= 0){
            stop("XiQ must be greater than zero.\n",call.=FALSE)
        }
    }else if(class(XiQ)=="matrix"){
        if(nrow(XiQ) != ((ncol(mydata)*p + 1)*ncol(mydata)) || ncol(XiQ) != ((ncol(mydata)*p + 1)*ncol(mydata))){
            stop("you have opted to give a full prior matrix for XiQ. However,\n", "it is not of the appropriate dimensions, which are ",(ncol(mydata)*p + 1)*ncol(mydata)," x ",(ncol(mydata)*p + 1)*ncol(mydata),".\n",call.=FALSE)
        }
    }else{
        stop("unrecognised form for XiQ.\n","XiQ must be a single numeric value or a matrix.\n",call.=FALSE)
    }
    #
    # Prior COV Degrees of Freedom
    if(class(gammaQ)=="NULL"){
        gammaQ <- tau + 1
    }else if(gammaQ <= 0){
        stop("gammaQ must be greater than zero.\n",call.=FALSE)
    }else if(gammaQ < tau){
        cat("Warning: you have selected a gammaQ less than tau, so gammaQ has been set to ",gammaQ + tau,".\n",sep="")
        gammaQ <- gammaQ + tau
    }else{gammaQ <- gammaQ}
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
    if(class(gammaS)=="NULL"){
        gammaS <- ncol(mydata) + 1
    }else if(gammaS <= ncol(mydata)){
        stop("with this data, the minimum value for gammaS is ",(ncol(mydata)+1),".\n",call.=FALSE)
    }else{gammaS <- gammaS}
    #
    #
    #
    return=list(coefprior=coefprior,XiBeta=XiBeta,XiQ=XiQ,XiSigma=XiSigma,gammaQ=gammaQ,gammaS=gammaS)
}

.bvartvpdata <- function(mydata,p,tau,timelab,irf.points){
    #
    Y <- as.matrix(mydata)
    #
    Tr <- as.numeric(nrow(Y))
    M <- as.numeric(ncol(Y))
    #
    K <- (M*p+1)*M
    #
    Ytau <- embed(Y,p+1); Ytau <- Ytau[,(M+1):ncol(Ytau)]; Ytau <- Ytau[(tau+1):(Tr-p),]
    #
    Z<-matrix(0,nrow=((Tr-tau-p)*M),ncol=K)
    for(i in 1:(Tr - tau - p)){
        Z[(((i-1)*M)+1):((i)*M),1:M] <- diag(M)
        for(j in 1:p){
            Z[(((i-1)*M)+1):((i)*M),(M+1+((M*M)*(j-1))):(M+(M*M)*(j))] <- t(kronecker(diag(M),Ytau[i,((j-1)*M+1):(j*M)]))
        }
    }
    #
    y <- t(Y[(tau+p+1):Tr,])
    #
    kT <- as.numeric(ncol(y))
    #
    if(class(timelab)=="NULL"){
        timelab <- (tau+p+1):Tr
    }else{
        timelab <- timelab[(tau+p+1):Tr]
    }
    #
    nIRFs <- 0
    if(class(irf.points) != "NULL"){
        nIRFs <- as.numeric(length(irf.points))
    }else{
        nIRFs <- Tr - (tau+p)
        irf.points <- timelab
    }
    #
    return=list(Y=Y,y=y,Z=Z,M=M,K=K,kT=kT,timelab=timelab,nIRFs=nIRFs,irf.points=irf.points)
}

.bvartvpprior <- function(Y,Z,p,M,K,kT,coefprior,tau,keep,irf.periods,XiBeta,XiQ,gammaQ,XiSigma,gammaS){
    #
    B0Pr <- matrix(NA,nrow=K,ncol=1); B0VPr <- 0; invB0VPr <- 0
    QPr <- 0; QVPr <- 0; QDraw <- 0; Qchol <- 0
    SPr <- 0; SVPr <- 0; SDraw <- 0; invSDraw <- 0
    #
    if(class(tau)!="NULL"){
        TauSamplingRun <- .Call("tsprior", Y,tau,M,K,p, PACKAGE = "BMR")
        #
        BetaVariance <- TauSamplingRun$BVPrOLS
        #
        B0Pr <- TauSamplingRun$BetaOLS
        B0VPr <- XiBeta*BetaVariance
        invB0VPr <- solve(B0VPr)
        #
        QPr <- XiQ*tau*BetaVariance
        QVPr <- gammaQ
        QDraw <- XiQ*diag(K)
        Qchol <- sqrt(XiQ)*diag(K)
        #
        if(class(XiSigma)=="numeric"){
            SPr <- diag(M)*XiSigma
        }else{SPr <- XiSigma}
        SVPr <- gammaS
        SDraw <- TauSamplingRun$SigmaOLS
        invSDraw <- solve(SDraw)
        #
    }else if(class(tau)=="NULL"){
        if(class(coefprior)=="numeric"){
            BPr<-rbind(rep(0,M),diag(M),matrix(rep(0,(p-1)*M*M),nrow=(p-1)*M,ncol=M))
            for(i in 1:M){BPr[(i+1),i]<-coefprior[i]}
        }else{BPr <- coefprior}
        #
        QSortIndex <- matrix(1,ncol=M,nrow=K/M)
        for(i in 1:M){
            QSortIndex[,i] <- QSortIndex[,i] + (i-1)
        }
        IndexCount <- M+1
        for(i in 1:p){
            for(k in 1:M){
                for(j in 1:M){
                    QSortIndex[M*(i-1)+1+j,k] <- IndexCount
                    IndexCount = IndexCount + 1
                }
            }
        }
        QSortIndex <- c(QSortIndex)
        #
        BPr <- c(t(BPr))[c(t(QSortIndex))]
        B0Pr <- matrix(BPr,ncol=1)
        #
        if(class(XiBeta)=="numeric"){
            B0VPr <- diag(K)*XiBeta
        }else{B0VPr <- XiBeta}
        invB0VPr <- solve(B0VPr)
        #
        if(class(XiQ)=="numeric"){
            QPr <- diag(K)*XiQ
        }else{QPr <- XiBeta}
        QVPr <- gammaQ
        QDraw <- QPr
        Qchol <- t(chol(QPr))
        #
        if(class(XiSigma)=="numeric"){
            SPr <- diag(M)*XiSigma
        }else{SPr <- XiSigma}
        SVPr <- gammaS
        SDraw <- SPr
        invSDraw <- solve(SDraw)
    }
    #
    return=list(B0Pr=B0Pr,B0VPr=B0VPr,invB0VPr=invB0VPr,QPr=QPr,QVPr=QVPr,SPr=SPr,SVPr=SVPr,QDraw=QDraw,Qchol=Qchol,SDraw=SDraw,invSDraw=invSDraw)
}

.bvartvpreps <- function(Y,y,Z,B0Pr,B0VPr,invB0VPr,QPr,QVPr,SPr,SVPr,QDraw,Qchol,SDraw,invSDraw,K,M,p,kT,keep,burnin,timelab,nIRFs,irf.points,irf.periods){
    #
    message('Starting Gibbs C++, ', date(),'.', sep="")
    RepsRun <- .Call("BVARTVPReps", y,Z,M,K,kT,keep,burnin,B0Pr,B0VPr,invB0VPr,QPr,QVPr,SPr,SVPr,QDraw,Qchol,SDraw,invSDraw, PACKAGE = "BMR")
    message('C++ reps finished, ', date(),'. Now generating IRFs.', sep="")
    #
    BetaDraws <- RepsRun$BetaDraws; QDraws <- RepsRun$QDraws; SDraws <- RepsRun$SDraws
    #
    # We need to reorder the matrices due to the stacking being not vectorised as we usually have it. Sigma should be okay.
    #
    QSortIndex <- matrix(1,ncol=M,nrow=K/M)
    for(i in 1:M){
        QSortIndex[,i] <- QSortIndex[,i] + (i-1)
    }
    IndexCount <- M+1
    for(i in 1:p){
        for(k in 1:M){
            for(j in 1:M){
                QSortIndex[M*(i-1)+1+j,k] <- IndexCount
                IndexCount = IndexCount + 1
            }
        }
    }
    QSortIndex <- c(QSortIndex)
    # BetaT is for use in IRFs later
    BetaT <- array(0,dim=c(K/M,M,keep,kT+1))
    #
    for(i in 1:keep){
        for(j in 1:(kT+1)){
            BetaDraws[,j,i] <- BetaDraws[QSortIndex,j,i]
            BetaT[,,i,j] <- matrix(BetaDraws[,j,i],ncol=M)
        }
        QDraws[,,i] <- QDraws[QSortIndex,QSortIndex,i]
    }
    #
    BetaMean <- apply(BetaDraws,c(1,2),mean)
    QMean <- apply(QDraws,c(1,2),mean)
    SigmaMean <- apply(SDraws,c(1,2),mean)
    #
    # IRF
    #
    Beta <- array(0,dim=c(K/M,M,keep))
    impmain <- array(0,dim=c(irf.periods,M,M,keep,nIRFs))
    KTC <- 1
    #
    for(i in 1:kT){
        if(KTC < (nIRFs+1)){
            if(timelab[i] == irf.points[KTC]){
                Beta <- BetaT[,,,i]
                ImpStore <- .Call("BVARTVPIRFs", M,K/M,keep,irf.periods,Beta,SDraws, PACKAGE = "BMR")
                ImpStore <- ImpStore$ImpStore
                ImpStoreE <- array(NA,dim=c(M,M,irf.periods,keep))
                for(j in 1:keep){
                    ImpStoreE[,,,j] <- ImpStore[,,((j-1)*irf.periods+1):(j*irf.periods)]
                }
                ImpStoreE<-aperm(ImpStoreE,c(3,1,2,4))
                #
                impmain[,,,,KTC] <- ImpStoreE
                #
                KTC <- KTC + 1
            }
        }
    }
    #
    impmain<-aperm(impmain,c(1,2,3,5,4))
    impfinal<-apply(impmain,c(1,2,3,4),sort)
    impfinal<-aperm(impfinal,c(2,3,4,5,1))
    #
    return=list(BetaMean=BetaMean,Betas=BetaT,QMean=QMean,QDraws=QDraws,SigmaMean=SigmaMean,SDraws=SDraws,IRFs=impfinal)
}