stationarity <- function(y,KPSSp=4,ADFp=8,print=TRUE){
  results <- .stationarity(y,KPSSp,ADFp,print)
  return=list(KPSS=results$KPSS, ADF=results$ADF, ADFLags=results$ADFLags)
}

.stationarity <- function(y,KPSSp=4,ADFp=8,print=TRUE){
  #
  n<-as.numeric(ncol(y))
  #
  DataNames <- as.character(c(colnames(y)))
  CritLevels <- as.character(c("1 Pct", "2.5 Pct", "5 Pct", "10 Pct"))
  ADFMain <- matrix(,nrow=3,ncol=n,dimnames=list(c("Time Trend:","Constant:","Neither:"),DataNames))
  #
  ADFModels <- as.character(c("Trend Model","Drift Model","None"))
  ADFLags <- matrix(,nrow=n,ncol=3,dimnames=list(DataNames,ADFModels))
  #
  for(i in 1:n){
    #ADF tests:
    ADF <- .KADF(as.matrix(y[,i],ncol=1),p=ADFp)
    #
    ress <- cbind(ADF$ADF)
    colnames(ress) <- DataNames[i]
    #
    ADFMain[,i] <- ress
    #
    numlags <- ADF$p
    ADFLags[i,] <- numlags
  }
  ADFMain <- data.frame(ADFMain)
  #
  ADFCrit <- .stationaritycrit(as.numeric(nrow(y)))
  colnames(ADFCrit) <- CritLevels
  #
  ADFMain <- cbind(ADFMain,ADFCrit)
  #
  numlags <- data.frame(ADFLags)
  colnames(numlags) <- c("Trend Model","Drift Model","None")
  #
  #KPSS Tests:
  KPSSMain<-matrix(,nrow=2,ncol=n,dimnames=list(c("Time Trend:","No Trend:"),DataNames))
  for(i in 1:n){
    KPSS <- .KKPSS(as.matrix(y[,i],ncol=1),p=KPSSp)
    #
    KPSSMain[,i] <- KPSS
  }
  #
  KPSSTrendCrit <- cbind(0.216,0.176,0.146,0.119)
  KPSSConsCrit <- cbind(0.739,0.574,0.463,0.347)
  KPSSCrit <- rbind(KPSSTrendCrit,KPSSConsCrit)
  #
  KPSSMain <- data.frame(KPSSMain)
  colnames(KPSSCrit) <- c("1 Pct","2.5 Pct","5 Pct","10 Pct")
  KPSSMain <- cbind(KPSSMain,KPSSCrit)
  #
  #
  if(print == T){
    cat("KPSS Tests:",KPSSp,"lags", "\n")
    print(KPSSMain)
    cat("", "\n")
    cat("Augmented Dickey-Fuller Tests:", "\n")
    print(ADFMain)
    cat("", "\n")
    cat("Number of Diff Lags for ADF Tests:", "\n")
    print(ADFLags)
  }
  return=list(KPSS=KPSSMain, ADF=ADFMain, ADFLags=ADFLags)
}

.KKPSS <- function(X,p=10){
  #
  Y<-as.matrix(X,ncol=1)
  #
  Epsilon <- matrix(NA,nrow=nrow(Y),ncol=2)
  X <- matrix(rep(1,nrow(Y))); X <- cbind(X,1:nrow(Y))
  for(i in 1:2){
    Alpha <- solve.qr(qr(as.matrix(X[,1:i],ncol=i)),Y)
    Epsilon[,i] <- Y - (as.matrix(X[,1:i],ncol=i))%*%Alpha
  }
  #
  KTest<-numeric(length=2)
  #
  for(i in 1:2){
    Residuals <- as.matrix(Epsilon[,i],ncol=1)
    kT <- nrow(Residuals)
    S = cumsum(Residuals); S = S^2; K = sum(S)/(kT^2)
    #
    KPSSfun <- function(p,kT,Residuals){
      S <- numeric(length=p)
      for(i in 1:p){
        ktemp <- t(matrix(Residuals[(1+i):nrow(Residuals),]))%*%matrix(Residuals[1:(nrow(Residuals)-i),])
        ktemp <- ktemp*1/kT
        S[i] <- (1 - (i/(1+p)))*ktemp
      }
      S <- 2*(sum(S))
      return(S)
    }
    #
    Sig <- (1/kT)*sum(Residuals^2) + KPSSfun(p,kT,Residuals)
    K <- K/Sig
    KTest[i] <- K
  }
  # Returns time trend result first, then constant
  return(KTest[2:1])
}

.stationaritycrit <- function(kT){
  #
  critarray <- array(NA,dim=c(6,4,3))
  # Tau 1
  critarray[,,3] <- cbind(c(-2.66, -2.62, -2.60, -2.58, -2.58, -2.58),
                          c(-2.26, -2.25, -2.24, -2.23, -2.23, -2.23),
                          c(-1.95, -1.95, -1.95, -1.95, -1.95, -1.95),
                          c(-1.60, -1.61, -1.61, -1.62, -1.62, -1.62))
  # Tau 2
  critarray[,,2] <- cbind(c(-3.75, -3.58, -3.51, -3.46, -3.44, -3.43),
                          c(-3.33, -3.22, -3.17, -3.14, -3.13, -3.12),
                          c(-3.00, -2.93, -2.89, -2.88, -2.87, -2.86),
                          c(-2.62, -2.60, -2.58, -2.57, -2.57, -2.57))
  # Tau 3
  critarray[,,1] <- cbind(c(-4.38, -4.15, -4.04, -3.99, -3.98, -3.96),
                          c(-3.95, -3.80, -3.73, -3.69, -3.68, -3.66),
                          c(-3.60, -3.50, -3.45, -3.43, -3.42, -3.41),
                          c(-3.24, -3.18, -3.15, -3.13, -3.13, -3.12))
  #
  critarray <- aperm(critarray,c(1,3,2))
  #returncrit <- matrix(NA,3,4)
  if(kT >= 500){returncrit <- critarray[6,,]}
  if(250 <= kT && kT < 500){returncrit <- critarray[5,,]}
  if(100 <= kT && kT < 250){returncrit <- critarray[4,,]}
  if(50 <= kT && kT < 100){returncrit <- critarray[3,,]}
  if(25 <= kT && kT < 50){returncrit <- critarray[2,,]}
  if(kT < 25){returncrit <- critarray[1,,]}
  #
  return(returncrit)
}

.KADF <- function(X,p=10){
  X <- as.matrix(X,ncol=1)
  XO <- matrix(X[-1,])
  #
  X <- embed(X,2) ; X <- X[,1] - X[,2]; X <- as.matrix(X,ncol=1)
  #
  XL <- embed(X,p+1); Y <- matrix(XL[,1]); XL <- XL[,2:ncol(XL)]
  kT <- as.numeric(nrow(Y))
  XL<-cbind(1:kT,rep(1,kT),XO[(p):(nrow(XO)-1),],XL)
  #
  whichp <- whichp <- .Call("ADFCheck", Y,XL,p,kT, PACKAGE = "BMR", DUP = FALSE)
  PMin <- whichp$BIC + 1
  #
  XL3<-embed(X,PMin[3]+1); Y3 <- matrix(XL3[,1]); XL3<-as.matrix(XL3[,2:ncol(XL3)],ncol=PMin[3]); XL3 <- cbind(XO[(PMin[3]):(nrow(XO)-1),],XL3)
  XL2<-embed(X,PMin[2]+1); Y2 <- matrix(XL2[,1]); XL2<-as.matrix(XL2[,2:ncol(XL2)],ncol=PMin[2]); XL2 <- cbind(rep(1,nrow(XL2)),XO[(PMin[2]):(nrow(XO)-1),],XL2)
  XL1<-embed(X,PMin[1]+1); Y1 <- matrix(XL1[,1]); XL1<-as.matrix(XL1[,2:ncol(XL1)],ncol=PMin[1]); XL1 <- cbind(rep(1,nrow(XL1)),1:nrow(XL1),XO[(PMin[1]):(nrow(XO)-1),],XL1)
  #
  TestStats<-numeric(length=3)
  #
  Beta <- solve.qr(qr(XL3),Y3); K <- length(Beta); Epsilon <- Y3 - XL3%*%Beta; 
  sigmaEST<-(1/(nrow(Y3)-K))*(t(Epsilon)%*%Epsilon); SE <- as.numeric(sigmaEST)*solve(t(XL3)%*%XL3); SE <- sqrt(diag(SE)); 
  TestStats[1] <- Beta[1,1]/SE[1]
  #
  Beta <- solve.qr(qr(XL2),Y2); K <- length(Beta); Epsilon <- Y2 - XL2%*%Beta; 
  sigmaEST<-(1/(nrow(Y2)-K))*(t(Epsilon)%*%Epsilon); SE <- as.numeric(sigmaEST)*solve(t(XL2)%*%XL2); SE <- sqrt(diag(SE)); 
  TestStats[2] <- Beta[2,1]/SE[2]
  #
  Beta <- solve.qr(qr(XL1),Y1); K <- length(Beta); Epsilon <- Y1 - XL1%*%Beta; 
  sigmaEST<-(1/(nrow(Y1)-K))*(t(Epsilon)%*%Epsilon); SE <- as.numeric(sigmaEST)*solve(t(XL1)%*%XL1); SE <- sqrt(diag(SE)); 
  TestStats[3] <- Beta[3,1]/SE[3]
  #
  return=list(ADF=TestStats[3:1],p=PMin)
}