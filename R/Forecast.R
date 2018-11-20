################################################################################
##
##   Copyright (C) 2011-2018 Keith O'Hara
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

forecast.Rcpp_bvarm <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

forecast.Rcpp_bvars <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

forecast.Rcpp_bvarcnw <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

forecast.Rcpp_bvarinw <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

forecast.Rcpp_cvar <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

forecast.Rcpp_dsge_gensys <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_dsge(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

forecast.Rcpp_dsge_uhlig <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_dsge(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

forecast.Rcpp_dsgevar_gensys <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

forecast.Rcpp_dsgevar_uhlig <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11,...)
{
    return=.forecast_var(obj,periods,shocks,plot,var_names,percentiles,use_mean,back_data,save,height,width)
}

#

.forecast_var <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11)
{    
    M <- obj$M
    n_draws <- dim(obj$beta_draws)[3]

    forecast_cube = obj$forecast(periods,shocks)$forecast_vals

    Y <- obj$Y

    #

    forecast_sorted <- apply(forecast_cube,c(1,2),sort)
    forecast_sorted <- aperm(forecast_sorted,c(2,3,1))
    forecast_mean <- apply(forecast_sorted,c(1,2),mean)
    
    upper_conf <- round(percentiles[3]*n_draws)
    mid_conf <- round(percentiles[2]*n_draws)
    lower_conf <- round(percentiles[1]*n_draws)

    #
    # Plot
    
    Time <- FCL <- FCU <- FCM <- NULL
    
    #
    
    plot_vals <- array(NA,dim=c((periods+back_data),4,M))
    
    for (i in 1:M) {
        FDataTemp <- 0
        
        if (use_mean == TRUE) { # Use the mean or middle percentile?
            FDataTemp <- data.frame(forecast_sorted[,i,lower_conf],forecast_mean[,i],forecast_sorted[,i,upper_conf])
        } else {
            FDataTemp <- data.frame(forecast_sorted[,i,lower_conf],forecast_sorted[,i,mid_conf],forecast_sorted[,i,upper_conf])
        }
        
        FDataTemp <- as.matrix(FDataTemp)
        
        #
        
        if(back_data > 0){
            FDataTemp <- rbind(matrix(rep(Y[(nrow(Y)-back_data+1):nrow(Y),i],3),ncol=3),FDataTemp)
        }
        
        FDataTemp <- cbind(FDataTemp,(as.numeric(nrow(Y))-back_data+1):(nrow(Y)+periods))
        plot_vals[,,i] <- FDataTemp
    }
    
    #

    if (plot)
    {
        vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
        
        if (save==TRUE) {
            if(class(dev.list()) != "NULL"){dev.off()}
            cairo_ps(filename="Forecast.eps",height=height,width=width)
        }
        
        #
        
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(M,1)))
        
        if (class(var_names) != "character") {
            var_names <- character(length=M)
            for (i in 1:M) {  
                var_names[i] <- paste("VAR",i,sep="")
            }
        }
        
        #
        
        if (back_data > 0) {
            # Include a dashed line to mark where the forecast begins
            for (i in 1:M) {
                FCastName <- var_names[i]
                FCDF <- plot_vals[,,i]
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
                FCDF <- plot_vals[,,i]
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

    return=list(forecast_mean=forecast_mean,plot_vals=plot_vals)
}

.forecast_dsge <- function(obj,periods=20,shocks=TRUE,plot=TRUE,var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,back_data=0,save=FALSE,height=13,width=11)
{
    forecast_cube = obj$forecast(periods,shocks)$forecast_vals

    M <- dim(forecast_cube)[2]
    n_draws <- dim(forecast_cube)[3]

    Y <- obj$estim_data

    #

    forecast_sorted <- apply(forecast_cube,c(1,2),sort)
    forecast_sorted <- aperm(forecast_sorted,c(2,3,1))
    forecast_mean <- apply(forecast_sorted,c(1,2),mean)
    
    upper_conf <- round(percentiles[3]*n_draws)
    mid_conf <- round(percentiles[2]*n_draws)
    lower_conf <- round(percentiles[1]*n_draws)

    #
    # Plot
    
    Time <- FCL <- FCU <- FCM <- NULL
    
    #
    
    plot_vals <- array(NA,dim=c((periods+back_data),4,M))
    
    for (i in 1:M) {
        FDataTemp <- 0
        
        if (use_mean == TRUE) { # Use the mean or middle percentile?
            FDataTemp <- data.frame(forecast_sorted[,i,lower_conf],forecast_mean[,i],forecast_sorted[,i,upper_conf])
        } else {
            FDataTemp <- data.frame(forecast_sorted[,i,lower_conf],forecast_sorted[,i,mid_conf],forecast_sorted[,i,upper_conf])
        }
        
        FDataTemp <- as.matrix(FDataTemp)
        
        #
        
        if(back_data > 0){
            FDataTemp <- rbind(matrix(rep(Y[(nrow(Y)-back_data+1):nrow(Y),i],3),ncol=3),FDataTemp)
        }
        
        FDataTemp <- cbind(FDataTemp,(as.numeric(nrow(Y))-back_data+1):(nrow(Y)+periods))
        plot_vals[,,i] <- FDataTemp
    }
    
    #

    if (plot)
    {
        vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
        
        if (save==TRUE) {
            if(class(dev.list()) != "NULL"){dev.off()}
            cairo_ps(filename="Forecast.eps",height=height,width=width)
        }
        
        #
        
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(M,1)))
        
        if (class(var_names) != "character") {
            var_names <- character(length=M)
            for (i in 1:M) {  
                var_names[i] <- paste("VAR",i,sep="")
            }
        }
        
        #
        
        if (back_data > 0) {
            # Include a dashed line to mark where the forecast begins
            for (i in 1:M) {
                FCastName <- var_names[i]
                FCDF <- plot_vals[,,i]
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
                FCDF <- plot_vals[,,i]
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
    
    return=list(forecast_mean=forecast_mean,plot_vals=plot_vals)
}
