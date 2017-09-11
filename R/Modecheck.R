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

modecheck.EDSGE <- function(obj,gridsize=200,scalepar=NULL,plottransform=FALSE,parnames=NULL,save=FALSE,height=13,width=13,...){
    .edsgemodecheck(obj,gridsize,scalepar,plottransform,parnames,save,height,width) 
}

modecheck.DSGEVAR <- function(obj,gridsize=200,scalepar=NULL,plottransform=FALSE,parnames=NULL,save=FALSE,height=13,width=13,...){
    .dsgevarmodecheck(obj,gridsize,scalepar,plottransform,parnames,save,height,width) 
}

.edsgemodecheck <- function(obj,gridsize=200,scalepar=NULL,plottransform=FALSE,parnames=NULL,save=FALSE,height=13,width=13){
    #
    dsgedata <- obj$data
    ObserveMat <- obj$ObserveMat
    partomats <- obj$partomats
    priorform <- obj$priorform
    priorpars <- obj$priorpars
    parbounds <- obj$parbounds
    #
    Parameters <- obj$Parameters
    nParam <- as.numeric(ncol(Parameters))
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
    priorform <- .edsgePrelimWork(dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)$priorform
    #
    parametersMode <- .DSGEParTransform(obj$parMode,priorform,parbounds,1)
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
    ParModeCheck <- array(NA,dim=c(gridsize+1,2,nParam))
    #
    ParamTemp <- parametersMode
    for(j in 1:nParam){
        ParamTemp <- parametersMode
        ParamGrid <- seq(from=(ParamTemp[j]-parametersModeHessian[j]),to=(ParamTemp[j]+parametersModeHessian[j]),length.out=(gridsize+1))
        for(i in 1:length(ParamGrid)){
            ParamTemp[j] <- ParamGrid[i]
            if(plottransform==FALSE){
                ParModeCheck[i,1,j] <- .DSGEParTransform(ParamTemp,priorform,parbounds,2)[j]
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
    plotpages <- 1
    if(plot==TRUE){
        ParameterVals <- LogPosterior <- NULL # CRAN check workaround
        #
        if(plottransform==FALSE){
            parametersMode <- .DSGEParTransform(parametersMode,priorform,parbounds,2)
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
        vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
        #
        ParamCount <- 1
        #
        for(ii in 1:plotpages){
            if(save==TRUE){
                if(class(dev.list()) != "NULL"){dev.off()}
                #
                SaveMode <- paste("DSGEModeCheck",as.character(ii),".eps",sep="")
                cairo_ps(filename=SaveMode,height=height,width=width)
            }
            grid.newpage()
            pushViewport(viewport(layout=grid.layout(MR,MC)))
            #
            for(i in 1:MR){
                for(k in 1:MC){
                    #
                    if(ParamCount <= (nParam)){
                        ParamDF <- data.frame(ParModeCheck[,,ParamCount])
                        colnames(ParamDF) <- c("ParameterVals","LogPosterior")
                        #
                        print(ggplot(data=(ParamDF),aes(x=ParameterVals)) + xlab("") + ylab("Log Posterior") + geom_vline(xintercept=parametersMode[ParamCount],linetype = "longdash") + geom_line(aes(y=LogPosterior),color="springgreen4") + labs(title=ParNames[ParamCount]) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
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

.dsgevarmodecheck <- function(obj,gridsize=200,scalepar=NULL,plottransform=TRUE,parnames=NULL,save=FALSE,height=13,width=13){
    #
    dsgedata <- obj$data
    ObserveMat <- obj$ObserveMat
    partomats <- obj$partomats
    priorform <- obj$priorform
    priorpars <- obj$priorpars
    parbounds <- obj$parbounds
    #
    Parameters <- obj$Parameters
    nParam <- as.numeric(ncol(Parameters))
    #
    lambda <- obj$lambda
    p <- obj$p
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
    prelimwork <- .dsgevarPrelimWork(dsgedata,lambda,p,obj$constant,matrix(0,1,1),FALSE,ObserveMat,partomats,priorform,priorpars,parbounds)
    kdata <- prelimwork$kdata; priorform <- prelimwork$priorform;
    dsgedata <- kdata$Y
    #
    parametersMode <- .DSGEParTransform(obj$parMode,priorform,parbounds,1)
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
    ParModeCheck <- array(NA,dim=c(gridsize+1,2,nParam))
    #
    ParamTemp <- parametersMode
    #
    if(is.finite(lambda)==TRUE){
        for(j in 1:nParam){
            ParamTemp <- parametersMode
            ParamGrid <- seq(from=(ParamTemp[j]-parametersModeHessian[j]),to=(ParamTemp[j]+parametersModeHessian[j]),length.out=(gridsize+1))
            for(i in 1:length(ParamGrid)){
                ParamTemp[j] <- ParamGrid[i]
                if(plottransform==FALSE){
                    ParModeCheck[i,1,j] <- .DSGEParTransform(ParamTemp,priorform,parbounds,2)[j]
                }else{
                    ParModeCheck[i,1,j] <- ParamGrid[i]
                }
                #
                ParModeCheck[i,2,j] <- (-1)*.DSGEVARLogPosterior(ParamTemp,kdata,lambda,p,kdata$YY,kdata$XY,kdata$XX,ObserveMat,partomats,priorform,priorpars,parbounds)
            }
        }
    }else{
        for(j in 1:nParam){
            ParamTemp <- parametersMode
            ParamGrid <- seq(from=(ParamTemp[j]-parametersModeHessian[j]),to=(ParamTemp[j]+parametersModeHessian[j]),length.out=(gridsize+1))
            for(i in 1:length(ParamGrid)){
                ParamTemp[j] <- ParamGrid[i]
                if(plottransform==FALSE){
                    ParModeCheck[i,1,j] <- .DSGEParTransform(ParamTemp,priorform,parbounds,2)[j]
                }else{
                    ParModeCheck[i,1,j] <- ParamGrid[i]
                }
                #
                ParModeCheck[i,2,j] <- (-1)*.DSGEVARLogPosteriorInf(ParamTemp,kdata,lambda,p,kdata$YY,kdata$XY,kdata$XX,ObserveMat,partomats,priorform,priorpars,parbounds)
            }
        }
    }
    #
    # Plot
    #
    plot=TRUE
    MR <- 0; MC <- 0
    plotpages <- 1
    if(plot==TRUE){
        ParameterVals <- LogPosterior <- NULL # CRAN check workaround
        #
        if(plottransform==FALSE){
            parametersMode <- .DSGEParTransform(parametersMode,priorform,parbounds,2)
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
        vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
        #
        ParamCount <- 1
        #
        for(ii in 1:plotpages){
            if(save==TRUE){
                if(class(dev.list()) != "NULL"){dev.off()}
                #
                SaveMode <- paste("DSGEVARModeCheck",as.character(ii),".eps",sep="")
                cairo_ps(filename=SaveMode,height=height,width=width)
            }
            grid.newpage()
            pushViewport(viewport(layout=grid.layout(MR,MC)))
            #
            for(i in 1:MR){
                for(k in 1:MC){
                    #
                    if(ParamCount <= (nParam)){
                        ParamDF <- data.frame(ParModeCheck[,,ParamCount])
                        colnames(ParamDF) <- c("ParameterVals","LogPosterior")
                        #
                        print(ggplot(data=(ParamDF),aes(x=ParameterVals)) + xlab("") + ylab("Log Posterior") + geom_vline(xintercept=parametersMode[ParamCount],linetype = "longdash") + geom_line(aes(y=LogPosterior),color="springgreen4") + labs(title=ParNames[ParamCount]) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
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
}