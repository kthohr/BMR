modecheck.EDSGE <- function(obj,scalepar=NULL,plottransform=TRUE,save=FALSE,height=13,width=13){
  .edsgemodecheck(obj,scalepar,plottransform,save,height,width) 
}

modecheck.DSGEVAR <- function(obj,scalepar=NULL,plottransform=TRUE,save=FALSE,height=13,width=13){
  .dsgevarmodecheck(obj,scalepar,plottransform,save,height,width) 
}

.edsgemodecheck <- function(obj,scalepar=NULL,plottransform=TRUE,save=FALSE,height=13,width=13){
  #
  dsgedata <- obj$data
  ObserveMat <- obj$ObserveMat
  partomats <- obj$partomats
  priorform <- obj$priorform
  priorpars <- obj$priorpars
  parbounds <- obj$parbounds
  #
  Parameters <- obj$Parameters
  ParNames <- colnames(Parameters)
  nParam <- as.numeric(ncol(Parameters))
  #
  parametersMode <- obj$ModeParamTrans
  parametersModeHessian <- obj$ModeHessian
  parametersModeHessian <- solve(parametersModeHessian)
  parametersModeHessian <- diag(parametersModeHessian)
  parametersModeHessian <- sqrt(parametersModeHessian)
  #
  if(class(scalepar)=="NULL"){
    scalepar <- obj$scalepar
  }
  parametersModeHessian <- parametersModeHessian*scalepar
  #
  GridSize <- 1000
  ParModeCheck <- array(NA,dim=c(GridSize+1,2,nParam))
  #
  ParamTemp <- parametersMode
  for(j in 1:nParam){
    ParamTemp <- parametersMode
    ParamGrid <- seq(from=(ParamTemp[j]-parametersModeHessian[j]),to=(ParamTemp[j]+parametersModeHessian[j]),length.out=(GridSize+1))
    for(i in 1:length(ParamGrid)){
      ParamTemp[j] <- ParamGrid[i]
      if(plottransform==FALSE){
        ParModeCheck[i,1,j] <- .DSGEParTransform(NULL,ParamTemp,priorform,parbounds)[j]
      }else{
        ParModeCheck[i,1,j] <- ParamGrid[i]
      }
      #
      ParModeCheck[i,2,j] <- (-1)*.dsgeposteriorfn(ParamTemp,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
    }
  }
  #
  # Plot
  #
  plot=TRUE
  MR <- 0; MC <- 0
  if(plot==TRUE){
    if(plottransform==FALSE){
      parametersMode <- .DSGEParTransform(NULL,parametersMode,priorform,parbounds)
    }
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
    }else if(nParam > 12 && nParam < 17){
      MR <- 4; MC <- 4
    }
    #
    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    #
    #
    if(save==TRUE){
      if(class(dev.list()) != "NULL"){dev.off()}
      cairo_ps(file="DSGEModeCheck.eps",height=height,width=width)
    }
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(MR,MC)))
    #
    ParamCount <- 1
    for(i in 1:MR){
      for(k in 1:MC){
        #
        if(ParamCount <= (nParam)){
          ParamDF <- data.frame(ParModeCheck[,,ParamCount])
          colnames(ParamDF) <- c("ParameterVals","LogPosterior")
          #
          print(ggplot(data=(ParamDF),aes(x=ParameterVals)) + xlab("Parameter Value") + ylab("Log Posterior") + geom_vline(xintercept=parametersMode[ParamCount],linetype = "longdash") + geom_line(aes(y=LogPosterior),color="springgreen4") + labs(title=ParNames[ParamCount]),vp = vplayout(i,k))
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
  #colnames(IRFs) <- varnames
  #
}

.dsgevarmodecheck <- function(obj,scalepar=NULL,plottransform=TRUE,save=FALSE,height=13,width=13){
  #
  require(ggplot2)
  require(grid)
  #
  dsgedata <- obj$data
  ObserveMat <- obj$ObserveMat
  partomats <- obj$partomats
  priorform <- obj$priorform
  priorpars <- obj$priorpars
  parbounds <- obj$parbounds
  #
  Parameters <- obj$Parameters
  ParNames <- colnames(Parameters)
  nParam <- as.numeric(ncol(Parameters))
  #
  parametersMode <- obj$ModeParamTrans
  parametersModeHessian <- obj$ModeHessian
  parametersModeHessian <- solve(parametersModeHessian)
  parametersModeHessian <- diag(parametersModeHessian)
  parametersModeHessian <- sqrt(parametersModeHessian)
  #
  lambda <- obj$lambda
  p <- (dim(obj$Beta)[1]/dim(obj$Beta)[2])
  #
  kdata <- .dsgevardata(dsgedata,p,FALSE)
  dsgedata <- kdata$Y
  #
  if(class(scalepar)=="NULL"){
    scalepar <- obj$scalepar
  }
  parametersModeHessian <- parametersModeHessian*scalepar
  #
  GridSize <- 1000
  ParModeCheck <- array(NA,dim=c(GridSize+1,2,nParam))
  #
  ParamTemp <- parametersMode
  #
  if(is.finite(lambda)==TRUE){
    for(j in 1:nParam){
      ParamTemp <- parametersMode
      ParamGrid <- seq(from=(ParamTemp[j]-parametersModeHessian[j]),to=(ParamTemp[j]+parametersModeHessian[j]),length.out=(GridSize+1))
      for(i in 1:length(ParamGrid)){
        ParamTemp[j] <- ParamGrid[i]
        if(plottransform==FALSE){
          ParModeCheck[i,1,j] <- .DSGEParTransform(NULL,ParamTemp,priorform,parbounds)[j]
        }else{
          ParModeCheck[i,1,j] <- ParamGrid[i]
        }
        #Remember, dsgeposteriorfn will return the NEGATIVE of the log posterior!
        ParModeCheck[i,2,j] <- (-1)*.DSGEVARLogPosterior(ParamTemp,kdata,lambda,p,kdata$YY,kdata$XY,kdata$XX,ObserveMat,partomats,priorform,priorpars,parbounds)
      }
    }
  }else{
    for(j in 1:nParam){
      ParamTemp <- parametersMode
      ParamGrid <- seq(from=(ParamTemp[j]-parametersModeHessian[j]),to=(ParamTemp[j]+parametersModeHessian[j]),length.out=(GridSize+1))
      for(i in 1:length(ParamGrid)){
        ParamTemp[j] <- ParamGrid[i]
        if(plottransform==FALSE){
          ParModeCheck[i,1,j] <- .DSGEParTransform(NULL,ParamTemp,priorform,parbounds)[j]
        }else{
          ParModeCheck[i,1,j] <- ParamGrid[i]
        }
        #Remember, dsgeposteriorfn will return the NEGATIVE of the log posterior!
        ParModeCheck[i,2,j] <- (-1)*.DSGEVARLogPosteriorInf(ParamTemp,kdata,lambda,p,kdata$YY,kdata$XY,kdata$XX,ObserveMat,partomats,priorform,priorpars,parbounds)
      }
    }
  }
  #
  # Plot
  #
  plot=TRUE
  MR <- 0; MC <- 0
  if(plot==TRUE){
    if(plottransform==FALSE){
      parametersMode <- .DSGEParTransform(NULL,parametersMode,priorform,parbounds)
    }
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
    }else if(nParam > 12 && nParam < 17){
      MR <- 4; MC <- 4
    }
    #
    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    #
    #
    if(save==TRUE){
      if(class(dev.list()) != "NULL"){dev.off()}
      cairo_ps(file="DSGEModeCheck.eps",height=height,width=width)
    }
    grid.newpage()
    pushViewport(viewport(layout=grid.layout(MR,MC)))
    #
    ParamCount <- 1
    for(i in 1:MR){
      for(k in 1:MC){
        #
        if(ParamCount <= (nParam)){
          ParamDF <- data.frame(ParModeCheck[,,ParamCount])
          colnames(ParamDF) <- c("ParameterVals","LogPosterior")
          #
          print(ggplot(data=(ParamDF),aes(x=ParameterVals)) + xlab("Parameter Value") + ylab("Log Posterior") + geom_vline(xintercept=parametersMode[ParamCount],linetype = "longdash") + geom_line(aes(y=LogPosterior),color="springgreen4") + labs(title=ParNames[ParamCount]),vp = vplayout(i,k))
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
  #colnames(IRFs) <- varnames
  #
}