plot.BVARM <- function(obj,save=FALSE,height=13,width=13){
  .plotbvarm(obj,save,height,width)
}

plot.BVARS <- function(obj,save=FALSE,height=13,width=13){
  .plotbvars(obj,save,height,width)
}

plot.BVARW <- function(obj,save=FALSE,height=13,width=13){
  .plotbvarw(obj,save,height,width)
}

plot.BVARTVP <- function(obj,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13){
  .plotbvartvp(obj,percentiles,save,height,width)
}

.plotbvarm <- function(obj,save=FALSE,height=13,width=13){
  Betas <- obj$BDraws
  #
  mydata <- obj$data
  constant <- obj$constant
  p <- floor(dim(Betas)[1]/dim(Betas)[2])
  M <- as.numeric(dim(Betas)[2])
  K <- as.numeric(dim(Betas)[1])
  keep<-as.numeric(dim(Betas)[3])
  # 
  BetaPerm <- aperm(Betas,c(3,1,2))
  #
  CoefLabels<-character(p*M)
  jj <- 1
  for(i in 1:p){
    for(k in 1:M){
      CoefLabels[jj]<-paste("L",i,".",k,sep="")
      jj = jj + 1
    }
  }
  #
  if(class(dev.list()) != "NULL"){dev.off()}
  #
  vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
  #
  CoefCount <- 1
  #
  if(constant==TRUE){
    for(i in 1:1){
      if(save==TRUE){cairo_ps(file="Constant.eps",height=(floor(height/M)),width=width)}
      pushViewport(viewport(layout=grid.layout(1,M)))
      #
      for(j in 1:M){
        VarName <- colnames(mydata)[j]
        #
        CFDF <- data.frame(BetaPerm[,i,j])
        colnames(CFDF) <- "CFDF"
        #
        if(j==1){
          print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Constant") + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
        }else{
          print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
        }
        Sys.sleep(0.6)
      }
      if(save==TRUE){dev.off()}
    }
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(file=SaveLag,height=height,width=width)}
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(M,M)))
      #
      for(j in 1:M){
        for(l in 1:M){
          VarName <- colnames(mydata)[l]
          #
          CFDF <- data.frame(BetaPerm[,((i-1)*M+j+1),l])
          colnames(CFDF) <- "CFDF"
          #
          #j==1 is for the variable title;l==1 is the coefficient label on y-axis
          if(j==1){
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }
          Sys.sleep(0.6)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }else{
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(file=SaveLag,height=height,width=width)}
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(M,M)))
      #
      for(j in 1:M){
        for(l in 1:M){
          VarName <- colnames(mydata)[l]
          #
          CFDF <- data.frame(BetaPerm[,((i-1)*M+j),l])
          colnames(CFDF) <- "CFDF"
          #
          #j==1 is for the variable title;l==1 is the coefficient label on y-axis
          if(j==1){
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }
          Sys.sleep(0.6)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
}

.plotbvars <- function(obj,save=FALSE,height=13,width=13){
  Psis <- obj$PDraws
  Betas <- obj$BDraws
  #
  mydata <- obj$data
  p <- floor(dim(Betas)[1]/dim(Betas)[2])
  M <- as.numeric(dim(Betas)[2])
  K <- as.numeric(dim(Betas)[1])
  keep<-as.numeric(dim(Betas)[3])
  # 
  BetaPerm <- aperm(Betas,c(3,1,2))
  PsiPerm <- aperm(Psis,c(3,1,2))
  #
  CoefLabels<-character(p*M)
  jj <- 1
  for(i in 1:p){
    for(k in 1:M){
      CoefLabels[jj]<-paste("L",i,".",k,sep="")
      jj = jj + 1
    }
  }
  #
  if(class(dev.list()) != "NULL"){dev.off()}
  #
  vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
  #
  CoefCount <- 1
  #
  for(i in 1:1){
    if(save==TRUE){cairo_ps(file="Psi.eps",height=(floor(height/M)),width=width)}
    pushViewport(viewport(layout=grid.layout(1,M)))
    #
    for(j in 1:M){
      VarName <- colnames(mydata)[j]
      #
      CFDF <- data.frame(PsiPerm[,1,j])
      colnames(CFDF) <- "CFDF"
      #
      if(j==1){
        print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Psi") + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
      }else{
        print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
      }
      Sys.sleep(0.6)
    }
    if(save==TRUE){dev.off()}
  }
  #
  for(i in 1:p){
    SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
    #
    if(save==TRUE){cairo_ps(file=SaveLag,height=height,width=width)}
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(M,M)))
    #
    for(j in 1:M){
      for(l in 1:M){
        VarName <- colnames(mydata)[l]
        #
        CFDF <- data.frame(BetaPerm[,((i-1)*M+j),l])
        colnames(CFDF) <- "CFDF"
        #
        #j==1 is for the variable title;l==1 is the coefficient label on y-axis
        if(j==1){
          if(l==1){
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            CoefCount <- CoefCount + 1
          }else{
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
          }
        }else{
          if(l==1){
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            CoefCount <- CoefCount + 1
          }else{
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
          }
        }
        Sys.sleep(0.6)
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
  plotSigma <- TRUE
  if(plotSigma==TRUE){
    Sigmas <- obj$SDraws
    SigmaPerm <- aperm(Sigmas,c(3,1,2))
    #
    SaveSig <- paste("Sigma.eps",sep="")
    if(save==TRUE){cairo_ps(file=SaveSig,height=height,width=width)}
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(M,M)))
    #
    for(j in 1:M){
      for(l in 1:M){
        VarNameY <- colnames(mydata)[j]
        VarNameX <- colnames(mydata)[l]
        #
        SDF <- data.frame(SigmaPerm[,j,l])
        colnames(SDF) <- "SDF"
        #
        #j==1 is for the variable title;l==1 is the coefficient label on y-axis
        if(j==1){
          if(l==1){
            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
          }else{
            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
          }
        }else{
          if(l==1){
            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
          }else{
            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
          }
        }
        Sys.sleep(0.6)
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
}

.plotbvarw <- function(obj,save=FALSE,height=13,width=13){
  Betas <- obj$BDraws
  #
  mydata <- obj$data
  constant <- obj$constant
  p <- floor(dim(Betas)[1]/dim(Betas)[2])
  M <- as.numeric(dim(Betas)[2])
  K <- as.numeric(dim(Betas)[1])
  keep <- as.numeric(dim(Betas)[3])
  # 
  BetaPerm <- aperm(Betas,c(3,1,2))
  #
  CoefLabels<-character(p*M)
  jj <- 1
  for(i in 1:p){
    for(k in 1:M){
      CoefLabels[jj]<-paste("L",i,".",k,sep="")
      jj = jj + 1
    }
  }
  #
  if(class(dev.list()) != "NULL"){dev.off()}
  #
  vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
  #
  CoefCount <- 1
  #
  if(constant==TRUE){
    for(i in 1:1){
      if(save==TRUE){cairo_ps(file="Constant.eps",height=(floor(height/M)),width=width)}
      pushViewport(viewport(layout=grid.layout(1,M)))
      #
      for(j in 1:M){
        VarName <- colnames(mydata)[j]
        #
        CFDF <- data.frame(BetaPerm[,i,j])
        colnames(CFDF) <- "CFDF"
        #
        if(j==1){
          print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Constant") + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
        }else{
          print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
        }
        Sys.sleep(0.6)
      }
      if(save==TRUE){dev.off()}
    }
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(file=SaveLag,height=height,width=width)}
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(M,M)))
      #
      for(j in 1:M){
        for(l in 1:M){
          VarName <- colnames(mydata)[l]
          #
          CFDF <- data.frame(BetaPerm[,((i-1)*M+j+1),l])
          colnames(CFDF) <- "CFDF"
          #
          #j==1 is for the variable title;l==1 is the coefficient label on y-axis
          if(j==1){
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }
          Sys.sleep(0.6)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }else{
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(file=SaveLag,height=height,width=width)}
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(M,M)))
      #
      for(j in 1:M){
        for(l in 1:M){
          VarName <- colnames(mydata)[l]
          #
          CFDF <- data.frame(BetaPerm[,((i-1)*M+j),l])
          colnames(CFDF) <- "CFDF"
          #
          #j==1 is for the variable title;l==1 is the coefficient label on y-axis
          if(j==1){
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }
          Sys.sleep(0.6)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
  plotSigma <- TRUE
  if(plotSigma==TRUE){
    Sigmas <- obj$SDraws
    SigmaPerm <- aperm(Sigmas,c(3,1,2))
    #
    SaveSig <- paste("Sigma.eps",sep="")
    if(save==TRUE){cairo_ps(file=SaveSig,height=height,width=width)}
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(M,M)))
    #
    for(j in 1:M){
      for(l in 1:M){
        VarNameY <- colnames(mydata)[j]
        VarNameX <- colnames(mydata)[l]
        #
        SDF <- data.frame(SigmaPerm[,j,l])
        colnames(SDF) <- "SDF"
        #
        #j==1 is for the variable title;l==1 is the coefficient label on y-axis
        if(j==1){
          if(l==1){
            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
          }else{
            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
          }
        }else{
          if(l==1){
            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
          }else{
            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
          }
        }
        Sys.sleep(0.6)
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
}

.plotbvartvp <- function(obj,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13){
  Betas <- obj$BDraws
  #
  mydata <- obj$data
  constant <- TRUE
  tau <- obj$tau
  #
  p <- floor(dim(Betas)[1]/dim(Betas)[2])
  M <- as.numeric(dim(Betas)[2])
  K <- as.numeric(dim(Betas)[1]) # K is (M*p)+1 in standard BVARs, but it's  ((M*p)+1)*M in BVARTVPs; let's stick with the standard BVAR definiton here.
  keep <- as.numeric(dim(Betas)[3])
  kT <- as.numeric(dim(Betas)[4])
  # kT x K x M x keep, for sorting
  BetaPerm <- aperm(Betas,c(4,1,2,3))
  BetaPerm2<-apply(BetaPerm,c(1,2,3),sort)
  BetaPerm <- BetaPerm2
  #
  CTPUpper <- round(percentiles[3]*keep)
  CTPMid <- round(percentiles[2]*keep)
  CTPLower <- round(percentiles[1]*keep)
  # 
  #
  CoefLabels<-character(p*M)
  jj <- 1
  for(i in 1:p){
    for(k in 1:M){
      CoefLabels[jj]<-paste("L",i,".",k,sep="")
      jj = jj + 1
    }
  }
  #
  if(class(dev.list()) != "NULL"){dev.off()}
  #
  vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
  #
  CoefCount <- 1
  #
  if(constant==TRUE){
    for(i in 1:1){
      if(save==TRUE){cairo_ps(file="Constant.eps",height=(floor(height/M)),width=width)}
      pushViewport(viewport(layout=grid.layout(1,M)))
      #
      for(j in 1:M){
        VarName <- colnames(mydata)[j]
        #
        CFDF <- data.frame(BetaPerm[CTPLower,,i,j],BetaPerm[CTPMid,,i,j],BetaPerm[CTPUpper,,i,j],tau:(kT+tau-1))
        colnames(CFDF) <- c("CTPL","CTPM","CTPU","Time")
        #
        if(j==1){
          print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab("Constant") + labs(title=paste(VarName)) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(1,j))
        }else{
          print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(NULL) + labs(title=paste(VarName)) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(1,j))
        }
        Sys.sleep(0.6)
      }
      if(save==TRUE){dev.off()}
    }
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(file=SaveLag,height=height,width=width)}
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(M,M)))
      #
      for(j in 1:M){
        for(l in 1:M){
          VarName <- colnames(mydata)[l]
          #
          CFDF <- data.frame(BetaPerm[CTPLower,,((i-1)*M+j+1),l],BetaPerm[CTPMid,,((i-1)*M+j+1),l],BetaPerm[CTPUpper,,((i-1)*M+j+1),l],tau:(kT+tau-1))
          colnames(CFDF) <- c("CTPL","CTPM","CTPU","Time")
          #CFDF <- data.frame(BetaPerm[,((i-1)*M+j+1),l])
          #colnames(CFDF) <- "CFDF"
          #
          #j==1 is for the variable title;l==1 is the coefficient label on y-axis
          if(j==1){
            if(l==1){
              print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1)+ geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(NULL) + labs(title=paste(VarName)) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(j,l))
            }
          }else{
            if(l==1){
              print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(paste(CoefLabels[CoefCount])) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(NULL) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(j,l))
            }
          }
          Sys.sleep(0.6)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }else{
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(file=SaveLag,height=height,width=width)}
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(M,M)))
      #
      for(j in 1:M){
        for(l in 1:M){
          VarName <- colnames(mydata)[l]
          #
          CFDF <- data.frame(BetaPerm[,((i-1)*M+j),l])
          colnames(CFDF) <- "CFDF"
          #
          #j==1 is for the variable title;l==1 is the coefficient label on y-axis
          if(j==1){
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(l==1){
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }
          Sys.sleep(0.6)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
}

plot.EDSGE <- function(obj,BinDenom=40,MCMCPlot=FALSE,save=FALSE,height=13,width=13){
  .plotedsge(obj,BinDenom,MCMCPlot,save,height,width)
}

plot.DSGEVAR <- function(obj,BinDenom=40,MCMCPlot=FALSE,save=FALSE,height=13,width=13){
  .plotdsgevar(obj,BinDenom,MCMCPlot,save,height,width)
}

.plotedsge <- function(obj,BinDenom=40,MCMCPlot=FALSE,save=FALSE,height=13,width=13){
  #
  Parameters <- obj$Parameters
  nParam <- as.numeric(ncol(Parameters))
  #
  if(nrow(Parameters)==0){
    stop("no MCMC draws detected.\n",call.=FALSE)
  }
  #
  if(class(dev.list()) != "NULL"){dev.off()}
  #
  vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
  #
  ParNames <- colnames(Parameters)
  if(class(ParNames) != "character"){
    ParNames <- character(length=nParam)
    for(i in 1:nParam){  
      ParNamess[i] <- paste("Parameter",i,sep="")
    }
  }
  #
  MR <- 0; MC <- 0
  if(nParam < 4){
    MR <- nParam; MC <- 1
  }else if(nParam == 4){
    MR <- 2; MC <-2
  }else if(nParam > 4 && nParam < 7){
    MR <- 3; MC <- 2
  }else if(nParam > 6 && nParam < 10){
    MR <- 3; MC <- 3
  }else if(nParam > 9 && nParam < 13){
    MR <- 4; MC <- 3
  }else if(nParam > 12 && nParam < 17){
    MR <- 4; MC <- 4
  }
  #
  nParGraphs<-1
  for(j in 1:nParGraphs){
    #
    if(save==TRUE){
      if(nParGraphs==1){
        cairo_ps(file="DSGEParameters.eps",height=height,width=width)
      }else{
        SaveParam <- paste("DSGEParameters_",j,".eps",sep="")
        #
        if(save==TRUE){cairo_ps(file=SaveParam,height=height,width=width)}
      }
    }
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(MR,MC)))
    #
    ParamCount <- 1
    for(i in 1:MR){
      for(k in 1:MC){
        #
        if(ParamCount <= nParam){
          ParamDF <- data.frame(Parameters[,ParamCount])
          colnames(ParamDF) <- c("DSGEPar")
          ParName <- ParNames[ParamCount]
          #
          ParamBin <- (max(ParamDF) - min(ParamDF))/BinDenom
          print(ggplot(data=(ParamDF),aes(DSGEPar)) + xlab("") + ylab("") + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=ParName),vp = vplayout(i,k))
          #print(ggplot(data=(ParamDF),aes(DSGEPar)) + xlab("") + ylab("") + geom_density(colour="darkblue",fill="purple",alpha=0.25) + labs(title=ParName),vp = vplayout(i,k))
          #
          ParamCount <- ParamCount + 1
          #
          Sys.sleep(0.6)
        }else{ParamCount <- ParamCount + 1}
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
  #
  #
  if(MCMCPlot==TRUE){
    for(j in 1:nParGraphs){
      #
      if(save==TRUE){
        if(class(dev.list()) != "NULL"){dev.off()}
        if(nParGraphs==1){
          cairo_ps(file="MCMCPlot.eps",height=height,width=width)
        }else{
          SaveParam <- paste("MCMCPlot_",j,".eps",sep="")
          #
          if(save==TRUE){cairo_ps(file=SaveParam,height=height,width=width)}
        }
      }
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(MR,MC)))
      #
      ParamCount <- 1
      for(i in 1:MR){
        for(k in 1:MC){
          #
          if(ParamCount <= nParam){
            ParamDF <- data.frame(Parameters[,ParamCount],1:nrow(Parameters))
            colnames(ParamDF) <- c("DSGEPar","Keep")
            ParName <- ParNames[ParamCount]
            #
            print(ggplot(data=(ParamDF),aes(x=Keep)) + xlab("Keep Run") + ylab(paste(ParName)) + geom_line(aes(y=DSGEPar),color="black"),vp = vplayout(i,k))
            #
            ParamCount <- ParamCount + 1
            #
            Sys.sleep(0.6)
          }else{ParamCount <- ParamCount + 1}
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
}

.plotdsgevar <- function(obj,BinDenom=40,MCMCPlot=FALSE,save=FALSE,height=13,width=13){
  #
  Parameters <- obj$Parameters
  nParam <- as.numeric(ncol(Parameters))
  #
  if(nrow(Parameters)==0){
    stop("no MCMC draws detected.\n",call.=FALSE)
  }
  #
  if(class(dev.list()) != "NULL"){dev.off()}
  #
  vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
  #
  ParNames <- colnames(Parameters)
  if(class(ParNames) != "character"){
    ParNames <- character(length=nParam)
    for(i in 1:nParam){  
      ParNamess[i] <- paste("Parameter",i,sep="")
    }
  }
  #
  MR <- 0; MC <- 0
  if(nParam < 4){
    MR <- nParam; MC <- 1
  }else if(nParam == 4){
    MR <- 2; MC <-2
  }else if(nParam > 4 && nParam < 7){
    MR <- 3; MC <- 2
  }else if(nParam > 6 && nParam < 10){
    MR <- 3; MC <- 3
  }else if(nParam > 9 && nParam < 13){
    MR <- 4; MC <- 3
  }else if(nParam > 12 && nParam < 17){
    MR <- 4; MC <- 4
  }
  #
  nParGraphs<-1
  for(j in 1:nParGraphs){
    #
    if(save==TRUE){
      if(nParGraphs==1){
        cairo_ps(file="DSGEParameters.eps",height=height,width=width)
      }else{
        SaveParam <- paste("DSGEParameters_",j,".eps",sep="")
        #
        if(save==TRUE){cairo_ps(file=SaveParam,height=height,width=width)}
      }
    }
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(MR,MC)))
    #
    ParamCount <- 1
    for(i in 1:MR){
      for(k in 1:MC){
        #
        if(ParamCount <= nParam){
          ParamDF <- data.frame(Parameters[,ParamCount])
          colnames(ParamDF) <- c("DSGEPar")
          ParName <- ParNames[ParamCount]
          #
          ParamBin <- (max(ParamDF) - min(ParamDF))/BinDenom
          print(ggplot(data=(ParamDF),aes(DSGEPar)) + xlab("") + ylab("") + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=ParName),vp = vplayout(i,k))
          #print(ggplot(data=(ParamDF),aes(DSGEPar)) + xlab("") + ylab("") + geom_density(colour="darkblue",fill="purple",alpha=0.25) + labs(title=ParName),vp = vplayout(i,k))
          #
          ParamCount <- ParamCount + 1
          #
          Sys.sleep(0.6)
        }else{ParamCount <- ParamCount + 1}
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
  #
  #
  if(MCMCPlot==TRUE){
    for(j in 1:nParGraphs){
      #
      if(save==TRUE){
        if(class(dev.list()) != "NULL"){dev.off()}
        if(nParGraphs==1){
          cairo_ps(file="MCMCPlot.eps",height=height,width=width)
        }else{
          SaveParam <- paste("MCMCPlot_",j,".eps",sep="")
          #
          if(save==TRUE){cairo_ps(file=SaveParam,height=height,width=width)}
        }
      }
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(MR,MC)))
      #
      ParamCount <- 1
      for(i in 1:MR){
        for(k in 1:MC){
          #
          if(ParamCount <= nParam){
            ParamDF <- data.frame(Parameters[,ParamCount],1:nrow(Parameters))
            colnames(ParamDF) <- c("DSGEPar","Keep")
            ParName <- ParNames[ParamCount]
            #
            print(ggplot(data=(ParamDF),aes(x=Keep)) + xlab("Keep Run") + ylab(paste(ParName)) + geom_line(aes(y=DSGEPar),color="black"),vp = vplayout(i,k))
            #
            ParamCount <- ParamCount + 1
            #
            Sys.sleep(0.6)
          }else{ParamCount <- ParamCount + 1}
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
}