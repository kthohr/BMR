################################################################################
##
##   Copyright (C) 2011-2017 Keith O'Hara
##
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

forecast.Rcpp_bvarm <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11,...){
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,useMean,backdata,save,height,width)
}

forecast.Rcpp_bvars <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11,...){
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,useMean,backdata,save,height,width)
}

forecast.Rcpp_bvarw <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11,...){
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,useMean,backdata,save,height,width)
}

forecast.Rcpp_cvar <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11,...){
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,useMean,backdata,save,height,width)
}

forecast.EDSGE <- function(obj,periods=20,plot=TRUE,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11,...){
    forecastdata <- .forecast_dsge(obj,periods,plot,percentiles,useMean,backdata,save,height,width)
    return=list(MeanForecast=forecastdata$MeanForecast,PointForecast=forecastdata$PointForecast,Forecasts=forecastdata$Forecasts)
}

forecast.DSGEVAR <- function(obj,periods=20,shocks=TRUE,plot=TRUE,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11,...){
    forecastdata <- .forecast_dsgevar(obj,periods,shocks,plot,percentiles,useMean,backdata,save,height,width)
    return=list(MeanForecast=forecastdata$MeanForecast,PointForecast=forecastdata$PointForecast,Forecasts=forecastdata$Forecasts)
}

.forecast_var <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11){
    
    #if(getRversion() >= "3.1.0") utils::suppressForeignCheck(names=c("Time", "FCL", "FCU", "FCM"))
    
    M <- obj$M
    n_draws <- dim(obj$beta_draws)[3]

    Y <- obj$Y

    forecast_cube = obj$forecast(periods,shocks)$forecast_vals

    #

    forecast_sorted <- apply(forecast_cube,c(1,2),sort)
    forecast_sorted <- aperm(forecast_sorted,c(2,3,1))
    forecast_mean <- apply(forecast_sorted,c(1,2),mean)
    
    upper_conf <- round(percentiles[3]*n_draws)
    mid_conf <- round(percentiles[2]*n_draws)
    lower_conf <- round(percentiles[1]*n_draws)

    #
    # Plotting
    #

    if (plot==TRUE) {
        Time <- FCL <- FCU <- FCM <- NULL

        #

        ForecastData <- array(NA,dim=c((periods+backdata),4,M))

        for (i in 1:M) {
            FDataTemp <- 0

            if (useMean == TRUE) { # Use the mean or middle percentile?
                FDataTemp <- data.frame(forecast_sorted[,i,lower_conf],forecast_mean[,i],forecast_sorted[,i,upper_conf])
            } else {
                FDataTemp <- data.frame(forecast_sorted[,i,lower_conf],forecast_sorted[,i,mid_conf],forecast_sorted[,i,upper_conf])
            }

            FDataTemp <- as.matrix(FDataTemp)
            
            #

            if(backdata > 0){
                FDataTemp <- rbind(matrix(rep(Y[(nrow(Y)-backdata+1):nrow(Y),i],3),ncol=3),FDataTemp)
            }
            
            FDataTemp <- cbind(FDataTemp,(as.numeric(nrow(Y))-backdata+1):(nrow(Y)+periods))
            ForecastData[,,i] <- FDataTemp
        }

        #

        vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}

        if(class(dev.list()) != "NULL"){dev.off()}

        if(save==TRUE){cairo_ps(filename="Forecast.eps",height=height,width=width)}

        pushViewport(viewport(layout=grid.layout(M,1)))

        if (class(var_names) != "character") {
            var_names <- character(length=M)
            for (i in 1:M) {  
                var_names[i] <- paste("VAR",i,sep="")
            }
        }
        
        #

        if (backdata > 0) {
            # Include a dashed line to mark where the forecast begins
            for (i in 1:M) {
                FCastName <- var_names[i]
                FCDF <- ForecastData[,,i]
                FCDF <- data.frame(FCDF)
                colnames(FCDF) <- c("FCL","FCM","FCU","Time")
                #
                print(ggplot(data=FCDF,aes(x=Time)) + xlab("Time") + ylab(paste("Forecast of ",FCastName)) + geom_ribbon(aes(ymin=FCL,ymax=FCU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0,colour='grey30') + geom_vline(xintercept=as.numeric(nrow(Y)),linetype = "longdash") + geom_line(aes(y=FCM),color="red",size=2) + theme(panel.background = element_rect(fill='white', colour='grey15')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,1))
                #
                Sys.sleep(0.6)
            }
        } else {
            for (i in 1:M) {
                FCastName <- var_names[i]
                FCDF <- ForecastData[,,i]
                FCDF <- data.frame(FCDF)
                colnames(FCDF) <- c("FCL","FCM","FCU","Time")
                #
                print(ggplot(data=FCDF,aes(x=Time)) + xlab("Time") + ylab(paste("Forecast of ",FCastName)) + geom_ribbon(aes(ymin=FCL,ymax=FCU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0,colour='grey30') + geom_line(aes(y=FCM),color="red",size=2) + theme(panel.background = element_rect(fill='white', colour='grey15')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,1))
                #
                Sys.sleep(0.6)
            }
        }
        if(save==TRUE){dev.off()}
    }
    #
    return=list(MeanForecast=forecast_mean,PointForecast=forecast_sorted[,,mid_conf])
}

.forecast_dsge <- function(obj,periods=20,plot=TRUE,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11){
    #
    DSGEPars <- obj$Parameters
    partomats <- obj$partomats
    ObserveMat <- obj$ObserveMat
    #
    Y <- obj$data
    Y <- as.matrix(Y,ncol=(ncol(Y)))
    M <- as.numeric(ncol(Y))
    p <- 1
    K <- 1 + p*M
    runs <- nrow(DSGEPars)
    #
    # Try to solve the model for the first set of values
    dsgemats <- partomats(DSGEPars[1,])
    dsgesolved <- SDSGE(dsgemats)
    StateMats <- statespace(dsgesolved)
    #
    IterMat <- diag(ncol(StateMats$F))
    Forecasts <- array(0,dim=c(periods,M,runs))
    #
    for(jj in 1:runs){
        dsgemats <- partomats(DSGEPars[jj,])
        dsgesolved <- SDSGE(dsgemats)
        StateMats <- statespace(dsgesolved)
        #
        IState <- .Call("DSGEKalmanFilt", Y,ObserveMat,dsgemats$ObsCons,StateMats$F,StateMats$G,dsgemats$shocks,dsgemats$MeasErrs,200, PACKAGE = "BMR")$dsgestate[,nrow(Y)]
        IState <- matrix(IState)
        #
        IterMat <- diag(ncol(StateMats$F))
        #
        for(kk in 1:periods){
            IterMat <- StateMats$F%*%IterMat
            Forecasts[kk,,jj] <- t(dsgemats$ObsCons + t(ObserveMat)%*%IterMat%*%IState)
        }
    }
    #
    ForecastsSorted <- apply(Forecasts,c(1,2),sort)
    ForecastsSorted <- aperm(ForecastsSorted,c(2,3,1))
    ForecastsMean <- apply(ForecastsSorted,c(1,2),mean)
    #
    UpperCInt <- round(percentiles[3]*runs)
    MidCInt <- round(percentiles[2]*runs)
    LowCInt <- round(percentiles[1]*runs)
    #
    # Plotting
    #
    if(plot==TRUE){
        Time <- FCL <- FCU <- FCM <- NULL # get around a CRAN check note
        #
        if(class(dev.list()) != "NULL"){dev.off()}
        #
        vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
        #
        ForecastData <- array(NA,dim=c((periods+backdata),4,M))
        for(i in 1:M){
            # Use the mean or selected percentile?
            if(useMean == T){
                FDataTemp<-data.frame(ForecastsSorted[,i,LowCInt],ForecastsMean[,i],ForecastsSorted[,i,UpperCInt])
            }else{
                FDataTemp<-data.frame(ForecastsSorted[,i,LowCInt],ForecastsSorted[,i,MidCInt],ForecastsSorted[,i,UpperCInt])
            }
            FDataTemp<-as.matrix(FDataTemp)
            #
            if(backdata > 0){
                FDataTemp<-rbind(matrix(rep(Y[(nrow(Y)-backdata+1):nrow(Y),i],3),ncol=3),FDataTemp)
            }
            FDataTemp<-cbind(FDataTemp,(as.numeric(nrow(Y))-backdata+1):(nrow(Y)+periods))
            ForecastData[,,i] <- FDataTemp
        }
        #
        if(class(dev.list()) != "NULL"){dev.off()}
        if(save==TRUE){cairo_ps(filename="Forecast.eps",height=height,width=width)}
        pushViewport(viewport(layout=grid.layout(M,1)))
        #
        # Include a dashed line to mark where the forecast begins
        if(backdata > 0){
            for(i in 1:M){
                FCastName <- colnames(obj$data)[i]
                FCDF <- ForecastData[,,i]
                FCDF <- data.frame(FCDF)
                colnames(FCDF) <- c("FCL","FCM","FCU","Time")
                #
                print(ggplot(data=FCDF,aes(x=Time)) + xlab("Time") + ylab(paste("Forecast of ",FCastName)) + geom_ribbon(aes(ymin=FCL,ymax=FCU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0,colour='grey30') + geom_vline(xintercept=as.numeric(nrow(Y)),linetype = "longdash") + geom_line(aes(y=FCM),color="red",size=2) + theme(panel.background = element_rect(fill='white', colour='grey15')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,1))
                #
                Sys.sleep(1)
            }
        }else{
            for(i in 1:M){
                FCastName <- colnames(obj$data)[i]
                FCDF <- ForecastData[,,i]
                FCDF <- data.frame(FCDF)
                colnames(FCDF) <- c("FCL","FCM","FCU","Time")
                #
                print(ggplot(data=FCDF,aes(x=Time)) + xlab("Time") + ylab(paste("Forecast of ",FCastName)) + geom_ribbon(aes(ymin=FCL,ymax=FCU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0,colour='grey30') + geom_line(aes(y=FCM),color="red",size=2) + theme(panel.background = element_rect(fill='white', colour='grey15')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,1))
                #
                Sys.sleep(1)
            }
        }
        if(save==TRUE){dev.off()}
    }
    #
    return=list(MeanForecast=ForecastsMean,PointForecast=ForecastsSorted[,,MidCInt],Forecasts=Forecasts)
}

.forecast_dsgevar <- function(obj,periods=20,shocks=TRUE,plot=TRUE,percentiles=c(.05,.50,.95),useMean=FALSE,backdata=0,save=FALSE,height=13,width=11){
    #
    Betas <- obj$Beta
    Sigmas <- obj$Sigma
    #
    if(class(Betas)=="NULL"){
        stop("no MCMC draws detected.\n",call.=FALSE)
    }
    #
    Y <- obj$data
    Y <- as.matrix(Y,ncol=(ncol(Y)))
    M <- as.numeric(ncol(Y))
    p <- obj$p
    K <- as.numeric(dim(Betas)[1])
    runs <- as.numeric(dim(Betas)[3])
    #
    constant <- obj$constant
    kcons <- 0
    if(constant==T){kcons <- 1}
    #
    inclshocks <- 0
    if(shocks==T){inclshocks <- 1}
    #
    kY0 <- embed(Y,p)
    kY0 <- kY0[nrow(kY0),]; kY0 <- as.matrix(kY0,ncol=1)
    if(constant==T){kY0 <- rbind(1,kY0)}
    kY0 <- t(kY0)
    #
    Forecasts <- .Call("dsgevarforecast", kY0,M,K,kcons,runs,periods,inclshocks,Betas,Sigmas, PACKAGE = "BMR")
    Forecasts <- Forecasts$Forecasts
    #
    ForecastsSorted <- apply(Forecasts,c(1,2),sort)
    ForecastsSorted <- aperm(ForecastsSorted,c(2,3,1))
    ForecastsMean <- apply(ForecastsSorted,c(1,2),mean)
    #
    UpperCInt <- round(percentiles[3]*runs)
    MidCInt <- round(percentiles[2]*runs)
    LowCInt <- round(percentiles[1]*runs)
    #
    # Plotting
    #
    if(plot==TRUE){
        Time <- FCL <- FCU <- FCM <- NULL # get around a CRAN check note
        #
        if(class(dev.list()) != "NULL"){dev.off()}
        #
        vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
        #
        ForecastData <- array(NA,dim=c((periods+backdata),4,M))
        for(i in 1:M){
            # Use the mean or selected percentile?
            if(useMean == T){
                FDataTemp<-data.frame(ForecastsSorted[,i,LowCInt],ForecastsMean[,i],ForecastsSorted[,i,UpperCInt])
            }else{
                FDataTemp<-data.frame(ForecastsSorted[,i,LowCInt],ForecastsSorted[,i,MidCInt],ForecastsSorted[,i,UpperCInt])
            }
            FDataTemp<-as.matrix(FDataTemp)
            #
            if(backdata > 0){
                FDataTemp<-rbind(matrix(rep(Y[(nrow(Y)-backdata+1):nrow(Y),i],3),ncol=3),FDataTemp)
            }
            FDataTemp<-cbind(FDataTemp,(as.numeric(nrow(Y))-backdata+1):(nrow(Y)+periods))
            ForecastData[,,i] <- FDataTemp
        }
        #
        if(class(dev.list()) != "NULL"){dev.off()}
        if(save==TRUE){cairo_ps(filename="Forecast.eps",height=height,width=width)}
        pushViewport(viewport(layout=grid.layout(M,1)))
        #
        # Include a dashed line to mark where the forecast begins
        if(backdata > 0){
            for(i in 1:M){
                FCastName <- colnames(obj$data)[i]
                FCDF <- ForecastData[,,i]
                FCDF <- data.frame(FCDF)
                colnames(FCDF) <- c("FCL","FCM","FCU","Time")
                #
                print(ggplot(data=FCDF,aes(x=Time)) + xlab("Time") + ylab(paste("Forecast of ",FCastName)) + geom_ribbon(aes(ymin=FCL,ymax=FCU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0,colour='grey30') + geom_vline(xintercept=as.numeric(nrow(Y)),linetype = "longdash") + geom_line(aes(y=FCM),color="red",size=2) + theme(panel.background = element_rect(fill='white', colour='grey15')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,1))
                #
                Sys.sleep(0.6)
            }
        }else{
            for(i in 1:M){
                FCastName <- colnames(obj$data)[i]
                FCDF <- ForecastData[,,i]
                FCDF <- data.frame(FCDF)
                colnames(FCDF) <- c("FCL","FCM","FCU","Time")
                #
                print(ggplot(data=FCDF,aes(x=Time)) + xlab("Time") + ylab(paste("Forecast of ",FCastName)) + geom_ribbon(aes(ymin=FCL,ymax=FCU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0,colour='grey30') + geom_line(aes(y=FCM),color="red",size=2) + theme(panel.background = element_rect(fill='white', colour='grey15')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,1))
                #
                Sys.sleep(0.6)
            }
        }
        if(save==TRUE){dev.off()}
    }
    #
    return=list(MeanForecast=ForecastsMean,PointForecast=ForecastsSorted[,,MidCInt],Forecasts=Forecasts)
}