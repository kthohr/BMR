states.EDSGE <- function(obj,percentiles=c(.05,.50,.95),varnames=NULL,useMean=FALSE,save=FALSE,height=13,width=11){
  #
  DSGEPars <- obj$Parameters
  partomats <- obj$partomats
  ObserveMat <- obj$ObserveMat
  #
  if(nrow(DSGEPars)==0){
    stop("no MCMC draws detected.\n",call.=FALSE)
  }
  #
  Y <- obj$data
  Y <- as.matrix(Y,ncol=(ncol(Y)))
  runs <- nrow(DSGEPars)
  #
  # Try to solve the model for the first set of values
  dsgemats <- partomats(DSGEPars[1,])
  dsgesolved <- SDSGE(dsgemats$A,dsgemats$B,dsgemats$C,dsgemats$D,dsgemats$F,dsgemats$G,dsgemats$H,dsgemats$J,dsgemats$K,dsgemats$L,dsgemats$M,dsgemats$N)
  StateMats <- .DSGEstatespace(dsgesolved$N,dsgesolved$P,dsgesolved$Q,dsgesolved$R,dsgesolved$S)
  nstates <- ncol(StateMats$F)
  #
  if(class(varnames) != "character"){
    varnames <- character(length=nstates)
    for(i in 1:nstates){  
      varnames[i] <- paste("State",i,sep="")
    }
  }
  #
  States <- array(0,dim=c(nrow(Y),nstates,runs))
  #
  for(jj in 1:runs){
    dsgemats <- partomats(DSGEPars[jj,])
    dsgesolved <- SDSGE(dsgemats$A,dsgemats$B,dsgemats$C,dsgemats$D,dsgemats$F,dsgemats$G,dsgemats$H,dsgemats$J,dsgemats$K,dsgemats$L,dsgemats$M,dsgemats$N)
    StateMats <- .DSGEstatespace(dsgesolved$N,dsgesolved$P,dsgesolved$Q,dsgesolved$R,dsgesolved$S)
    #
    States[,,jj] <- t(.Call("DSGEKalmanFilt", Y,ObserveMat,dsgemats$ObsCons,StateMats$F,StateMats$G,dsgesolved$N,dsgemats$shocks,dsgemats$MeasErrs,200, PACKAGE = "BMR", DUP = FALSE)$dsgestate)
  }
  #
  StatesSorted <- apply(States,c(1,2),sort)
  StatesSorted <- aperm(StatesSorted,c(2,3,1))
  StatesMean <- apply(StatesSorted,c(1,2),mean)
  #
  UpperCInt <- round(percentiles[3]*runs)
  MidCInt <- round(percentiles[2]*runs)
  LowCInt <- round(percentiles[1]*runs)
  #
  StateData <- array(NA,dim=c(nrow(Y),4,nstates))
  for(i in 1:nstates){
    # Use the mean or selected percentile?
    if(useMean == T){
      FDataTemp<-data.frame(StatesSorted[,i,LowCInt],StatesMean[,i],StatesSorted[,i,UpperCInt])
    }else{
      FDataTemp<-data.frame(StatesSorted[,i,LowCInt],StatesSorted[,i,MidCInt],StatesSorted[,i,UpperCInt])
    }
    FDataTemp <- as.matrix(FDataTemp)
    FDataTemp <- cbind(FDataTemp,1:nrow(Y))
    #
    StateData[,,i] <- FDataTemp
  }
  #
  #
  #
  if(class(dev.list()) != "NULL"){dev.off()}
  #
  vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
  #
  MR <- 0; MC <- 0
  plotpages <- 1
  #
  if(nstates < 4){
    MR <- nstates; MC <- 1
  }else if(nstates == 4){
    MR <- 2; MC <-2
  }else if(nstates > 4 && nstates < 7){
    MR <- 3; MC <- 2
  }else if(nstates > 6 && nstates < 10){
    MR <- 3; MC <- 3
  }else if(nstates > 9 && nstates < 13){
    MR <- 4; MC <- 3
  }else if(nstates > 12 && nstates < 25){
    MR <- 4; MC <- 3
    plotpages <- 2
  }else if(nstates > 24 && nstates < 37){
    MR <- 4; MC <- 3
    plotpages <- 3
  }else if(nstates > 36 && nstates < 49){
    MR <- 4; MC <- 3
    plotpages <- 4
  }else if(nstates > 48 && nstates < 61){
    MR <- 4; MC <- 3
    plotpages <- 5
  }else if(nstates > 60 && nstates < 73){
    MR <- 4; MC <- 3
    plotpages <- 6
  }
  #
  StateCount <- 1
  for(j in 1:plotpages){
    #
    if(save==TRUE){
      if(plotpages==1){
        cairo_ps(file="States.eps",height=height,width=width)
      }else{
        SaveState <- paste("States_",j,".eps",sep="")
        #
        if(save==TRUE){cairo_ps(file=SaveState,height=height,width=width)}
      }
    }
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(MR,MC)))
    #
    for(i in 1:MR){
      for(k in 1:MC){
        #
        if(StateCount <= nstates){
          StateName <- varnames[StateCount]
          SDF <- StateData[,,StateCount]
          SDF <- data.frame(SDF)
          colnames(SDF) <- c("SL","SM","SU","Time")
          #
          print(ggplot(data=SDF,aes(x=Time)) + xlab("Time") + ylab(paste(StateName)) + geom_ribbon(aes(ymin=SL,ymax=SU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_line(aes(y=SM),color="red",size=1) + theme(panel.background = element_rect(fill='white', colour='grey15')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
          #
          StateCount <- StateCount + 1
          #
          Sys.sleep(0.6)
          #
        }else{StateCount <- StateCount + 1}
      }
    }
    if(save==TRUE){dev.off()}
    #
  }
  #
  return=list(MeanState=StatesMean,States=StateData)
}