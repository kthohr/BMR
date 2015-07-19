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

plot.BVARM <- function(x,type=1,save=FALSE,height=13,width=13,...){
  .plotbvarm(x,type,save,height,width)
}

plot.BVARS <- function(x,type=1,plotSigma=TRUE,save=FALSE,height=13,width=13,...){
  .plotbvars(x,type,plotSigma,save,height,width)
}

plot.BVARW <- function(x,type=1,plotSigma=TRUE,save=FALSE,height=13,width=13,...){
  .plotbvarw(x,type,plotSigma,save,height,width)
}

plot.BVARTVP <- function(x,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13,...){
  .plotbvartvp(x,percentiles,save,height,width)
}

plot.EDSGE <- function(x,parnames=NULL,BinDenom=40,MCMCplot=FALSE,save=FALSE,height=13,width=13,...){
  .plotedsge(x,parnames,BinDenom,MCMCplot,save,height,width)
}

plot.DSGEVAR <- function(x,parnames=NULL,BinDenom=40,MCMCplot=FALSE,save=FALSE,height=13,width=13,...){
  .plotdsgevar(x,parnames,BinDenom,MCMCplot,save,height,width)
}

.plotbvarm <- function(obj,type=1,save=FALSE,height=13,width=13){
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
  BinDenom <- 40
  ParamBin <- 1
  #
  if(constant==TRUE){
    for(i in 1:1){
      if(save==TRUE){cairo_ps(filename="Constant.eps",height=(floor(height/M)),width=width)}
      pushViewport(viewport(layout=grid.layout(1,M)))
      #
      for(j in 1:M){
        VarName <- colnames(mydata)[j]
        #
        CFDF <- data.frame(BetaPerm[,i,j])
        colnames(CFDF) <- "CFDF"
        #
        if(j==1){
          if(type==1){
            ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Constant") + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
          }else{
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Constant") + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
          }
        }else{
          if(type==1){
            ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
          }else{
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
          }
        }
        Sys.sleep(0.3)
      }
      if(save==TRUE){dev.off()}
    }
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(filename=SaveLag,height=height,width=width)}
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
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }
            }else{
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              }
            }
          }else{
            if(l==1){
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }
            }else{
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              }
            }
          }
          Sys.sleep(0.3)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }else{
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(filename=SaveLag,height=height,width=width)}
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
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }
            }else{
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              }
            }
          }else{
            if(l==1){
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }
            }else{
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              }
            }
          }
          Sys.sleep(0.3)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
}

.plotbvars <- function(obj,type=1,plotSigma=TRUE,save=FALSE,height=13,width=13){
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
  BinDenom <- 40
  ParamBin <- 1
  #
  for(i in 1:1){
    if(save==TRUE){cairo_ps(filename="Psi.eps",height=(floor(height/M)),width=width)}
    pushViewport(viewport(layout=grid.layout(1,M)))
    #
    for(j in 1:M){
      VarName <- colnames(mydata)[j]
      #
      CFDF <- data.frame(PsiPerm[,1,j])
      colnames(CFDF) <- "CFDF"
      #
      if(j==1){
        if(type==1){
          ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
          print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Psi") + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
        }else{
          print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Psi") + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
        }
      }else{
        if(type==1){
          ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
          print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
        }else{
          print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
        }
      }
      Sys.sleep(0.3)
    }
    if(save==TRUE){dev.off()}
  }
  #
  for(i in 1:p){
    SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
    #
    if(save==TRUE){cairo_ps(filename=SaveLag,height=height,width=width)}
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
            if(type==1){
              ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }
          }else{
            if(type==1){
              ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }
        }else{
          if(l==1){
            if(type==1){
              ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              CoefCount <- CoefCount + 1
            }
          }else{
            if(type==1){
              ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
            }
          }
        }
        Sys.sleep(0.3)
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
  if(plotSigma==TRUE){
    Sigmas <- obj$SDraws
    SigmaPerm <- aperm(Sigmas,c(3,1,2))
    #
    SaveSig <- paste("Sigma.eps",sep="")
    if(save==TRUE){cairo_ps(filename=SaveSig,height=height,width=width)}
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
            if(type==1){
              ParamBin <- (max(SDF) - min(SDF))/BinDenom
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + labs(title=paste(VarNameX)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(type==1){
              ParamBin <- (max(SDF) - min(SDF))/BinDenom
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + labs(title=paste(VarNameX)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
            }
          }
        }else{
          if(l==1){
            if(type==1){
              ParamBin <- (max(SDF) - min(SDF))/BinDenom
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(type==1){
              ParamBin <- (max(SDF) - min(SDF))/BinDenom
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
            }
          }
        }
        Sys.sleep(0.3)
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
}

.plotbvarw <- function(obj,type=1,plotSigma=TRUE,save=FALSE,height=13,width=13){
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
  BinDenom <- 40
  ParamBin <- 1
  #
  if(constant==TRUE){
    for(i in 1:1){
      if(save==TRUE){cairo_ps(filename="Constant.eps",height=(floor(height/M)),width=width)}
      pushViewport(viewport(layout=grid.layout(1,M)))
      #
      for(j in 1:M){
        VarName <- colnames(mydata)[j]
        #
        CFDF <- data.frame(BetaPerm[,i,j])
        colnames(CFDF) <- "CFDF"
        #
        if(j==1){
          if(type==1){
            ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Constant") + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
          }else{
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Constant") + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
          }
        }else{
          if(type==1){
            ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
          }else{
            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
          }
        }
        Sys.sleep(0.3)
      }
      if(save==TRUE){dev.off()}
    }
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(filename=SaveLag,height=height,width=width)}
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
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }
            }else{
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              }
            }
          }else{
            if(l==1){
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }
            }else{
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              }
            }
          }
          Sys.sleep(0.3)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }else{
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(filename=SaveLag,height=height,width=width)}
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
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }
            }else{
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              }
            }
          }else{
            if(l==1){
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                CoefCount <- CoefCount + 1
              }
            }else{
              if(type==1){
                ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
              }else{
                print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
              }
            }
          }
          Sys.sleep(0.3)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
  if(plotSigma==TRUE){
    Sigmas <- obj$SDraws
    SigmaPerm <- aperm(Sigmas,c(3,1,2))
    #
    SaveSig <- paste("Sigma.eps",sep="")
    if(save==TRUE){cairo_ps(filename=SaveSig,height=height,width=width)}
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
            if(type==1){
              ParamBin <- (max(SDF) - min(SDF))/BinDenom
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + labs(title=paste(VarNameX)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(type==1){
              ParamBin <- (max(SDF) - min(SDF))/BinDenom
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + labs(title=paste(VarNameX)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
            }
          }
        }else{
          if(l==1){
            if(type==1){
              ParamBin <- (max(SDF) - min(SDF))/BinDenom
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
            }
          }else{
            if(type==1){
              ParamBin <- (max(SDF) - min(SDF))/BinDenom
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
            }else{
              print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
            }
          }
        }
        Sys.sleep(0.3)
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
  K <- as.numeric(dim(Betas)[1])
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
  Time <- CTPL <- CTPM <- CTPU <- NULL
  #
  CoefCount <- 1
  #
  if(constant==TRUE){
    for(i in 1:1){
      if(save==TRUE){cairo_ps(filename="Constant.eps",height=(floor(height/M)),width=width)}
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
        Sys.sleep(0.3)
      }
      if(save==TRUE){dev.off()}
    }
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(filename=SaveLag,height=height,width=width)}
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(M,M)))
      #
      for(j in 1:M){
        for(l in 1:M){
          VarName <- colnames(mydata)[l]
          #
          CFDF <- data.frame(BetaPerm[CTPLower,,((i-1)*M+j+1),l],BetaPerm[CTPMid,,((i-1)*M+j+1),l],BetaPerm[CTPUpper,,((i-1)*M+j+1),l],tau:(kT+tau-1))
          colnames(CFDF) <- c("CTPL","CTPM","CTPU","Time")
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
          Sys.sleep(0.3)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }else{
    for(i in 1:p){
      SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
      #
      if(save==TRUE){cairo_ps(filename=SaveLag,height=height,width=width)}
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
          Sys.sleep(0.3)
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
}

.plotedsge <- function(obj,parnames=NULL,BinDenom=40,MCMCplot=FALSE,save=FALSE,height=13,width=13){
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
  ParNames <- parnames
  if(is.null(parnames)==TRUE){
    ParNames <- colnames(Parameters)
    if(class(ParNames) != "character"){
      ParNames <- character(length=nParam)
      for(i in 1:nParam){  
        ParNames[i] <- paste("Parameter",i,sep="")
      }
    }
  }
  #
  MR <- 0; MC <- 0
  plotpages <- 1
  #
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
  }else if(nParam > 12 && nParam < 25){
    MR <- 4; MC <- 3
    plotpages <- 2
  }else if(nParam > 24 && nParam < 37){
    MR <- 4; MC <- 3
    plotpages <- 3
  }else if(nParam > 36 && nParam < 49){
    MR <- 4; MC <- 3
    plotpages <- 4
  }else if(nParam > 48 && nParam < 61){
    MR <- 4; MC <- 3
    plotpages <- 5
  }else if(nParam > 60 && nParam < 73){
    MR <- 4; MC <- 3
    plotpages <- 6
  }
  #
  ParamCount <- 1
  for(j in 1:plotpages){
    #
    if(save==TRUE){
      if(plotpages==1){
        cairo_ps(filename="DSGEParameters.eps",height=height,width=width)
      }else{
        SaveParam <- paste("DSGEParameters_",j,".eps",sep="")
        #
        if(save==TRUE){cairo_ps(filename=SaveParam,height=height,width=width)}
      }
    }
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(MR,MC)))
    #
    DSGEPar <- Keep <- NULL # CRAN check workaround
    #
    for(i in 1:MR){
      for(k in 1:MC){
        #
        if(ParamCount <= nParam){
          ParamDF <- data.frame(Parameters[,ParamCount])
          colnames(ParamDF) <- c("DSGEPar")
          ParName <- ParNames[ParamCount]
          #
          ParamBin <- (max(ParamDF) - min(ParamDF))/BinDenom
          print(ggplot(data=(ParamDF),aes(DSGEPar)) + xlab("") + ylab("") + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=ParName) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
          #
          ParamCount <- ParamCount + 1
          #
          Sys.sleep(0.3)
        }else{ParamCount <- ParamCount + 1}
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
  #
  #
  if(MCMCplot==TRUE){
    ParamCount <- 1
    for(j in 1:plotpages){
      #
      if(save==TRUE){
        if(class(dev.list()) != "NULL"){dev.off()}
        if(plotpages==1){
          cairo_ps(filename="MCMCplot.eps",height=height,width=width)
        }else{
          SaveParam <- paste("MCMCplot_",j,".eps",sep="")
          #
          if(save==TRUE){cairo_ps(filename=SaveParam,height=height,width=width)}
        }
      }
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(MR,MC)))
      #
      for(i in 1:MR){
        for(k in 1:MC){
          #
          if(ParamCount <= nParam){
            ParamDF <- data.frame(Parameters[,ParamCount],1:nrow(Parameters))
            colnames(ParamDF) <- c("DSGEPar","Keep")
            ParName <- ParNames[ParamCount]
            #
            print(ggplot(data=(ParamDF),aes(x=Keep)) + xlab("Keep Run") + ylab(paste(ParName)) + geom_line(aes(y=DSGEPar),color="black") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
            #
            ParamCount <- ParamCount + 1
            #
            Sys.sleep(0.3)
          }else{ParamCount <- ParamCount + 1}
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
}

.plotdsgevar <- function(obj,parnames=NULL,BinDenom=40,MCMCplot=FALSE,save=FALSE,height=13,width=13){
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
  ParNames <- parnames
  if(is.null(parnames)==TRUE){
    ParNames <- colnames(Parameters)
    if(class(ParNames) != "character"){
      ParNames <- character(length=nParam)
      for(i in 1:nParam){  
        ParNames[i] <- paste("Parameter",i,sep="")
      }
    }
  }
  #
  MR <- 0; MC <- 0
  plotpages <- 1
  #
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
  }else if(nParam > 12 && nParam < 25){
    MR <- 4; MC <- 3
    plotpages <- 2
  }else if(nParam > 24 && nParam < 37){
    MR <- 4; MC <- 3
    plotpages <- 3
  }else if(nParam > 36 && nParam < 49){
    MR <- 4; MC <- 3
    plotpages <- 4
  }else if(nParam > 48 && nParam < 61){
    MR <- 4; MC <- 3
    plotpages <- 5
  }else if(nParam > 60 && nParam < 73){
    MR <- 4; MC <- 3
    plotpages <- 6
  }
  #
  ParamCount <- 1
  for(j in 1:plotpages){
    #
    if(save==TRUE){
      if(plotpages==1){
        cairo_ps(filename="DSGEParameters.eps",height=height,width=width)
      }else{
        SaveParam <- paste("DSGEParameters_",j,".eps",sep="")
        #
        if(save==TRUE){cairo_ps(filename=SaveParam,height=height,width=width)}
      }
    }
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(MR,MC)))
    #
    DSGEPar <- Keep <- NULL # CRAN check workaround
    #
    for(i in 1:MR){
      for(k in 1:MC){
        #
        if(ParamCount <= nParam){
          ParamDF <- data.frame(Parameters[,ParamCount])
          colnames(ParamDF) <- c("DSGEPar")
          ParName <- ParNames[ParamCount]
          #
          ParamBin <- (max(ParamDF) - min(ParamDF))/BinDenom
          print(ggplot(data=(ParamDF),aes(DSGEPar)) + xlab("") + ylab("") + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=ParName) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
          #
          ParamCount <- ParamCount + 1
          #
          Sys.sleep(0.3)
        }else{ParamCount <- ParamCount + 1}
      }
    }
    if(save==TRUE){dev.off()}
  }
  #
  #
  #
  if(MCMCplot==TRUE){
    ParamCount <- 1
    for(j in 1:plotpages){
      #
      if(save==TRUE){
        if(class(dev.list()) != "NULL"){dev.off()}
        if(plotpages==1){
          cairo_ps(filename="MCMCplot.eps",height=height,width=width)
        }else{
          SaveParam <- paste("MCMCplot_",j,".eps",sep="")
          #
          if(save==TRUE){cairo_ps(filename=SaveParam,height=height,width=width)}
        }
      }
      grid.newpage()
      pushViewport(viewport(layout=grid.layout(MR,MC)))
      #
      for(i in 1:MR){
        for(k in 1:MC){
          #
          if(ParamCount <= nParam){
            ParamDF <- data.frame(Parameters[,ParamCount],1:nrow(Parameters))
            colnames(ParamDF) <- c("DSGEPar","Keep")
            ParName <- ParNames[ParamCount]
            #
            print(ggplot(data=(ParamDF),aes(x=Keep)) + xlab("Keep Run") + ylab(paste(ParName)) + geom_line(aes(y=DSGEPar),color="black") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
            #
            ParamCount <- ParamCount + 1
            #
            Sys.sleep(0.3)
          }else{ParamCount <- ParamCount + 1}
        }
      }
      if(save==TRUE){dev.off()}
    }
  }
  #
}