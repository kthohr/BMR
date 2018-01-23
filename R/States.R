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

states.Rcpp_dsge_gensys <- function(obj,percentiles=c(.05,.50,.95),var_names=NULL,use_mean=FALSE,save=FALSE,height=13,width=11,...)
{
    return=.states_int(obj,percentiles,var_names,use_mean,save,height,width)
}

states.Rcpp_dsge_uhlig <- function(obj,percentiles=c(.05,.50,.95),var_names=NULL,use_mean=FALSE,save=FALSE,height=13,width=11,...)
{
    return=.states_int(obj,percentiles,var_names,use_mean,save,height,width)
}

states.Rcpp_dsgevar_gensys <- function(obj,percentiles=c(.05,.50,.95),var_names=NULL,use_mean=FALSE,save=FALSE,height=13,width=11,...)
{
    return=.states_int(obj,percentiles,var_names,use_mean,save,height,width)
}

states.Rcpp_dsgevar_uhlig <- function(obj,percentiles=c(.05,.50,.95),var_names=NULL,use_mean=FALSE,save=FALSE,height=13,width=11,...)
{
    return=.states_int(obj,percentiles,var_names,use_mean,save,height,width)
}

#

.states_int <- function(obj,percentiles=c(.05,.50,.95),var_names=NULL,use_mean=FALSE,save=FALSE,height=13,width=11,...)
{
    
    n_draws <- dim(obj$dsge_draws)[1]

    if (n_draws <= 0) {
        stop("error: no MCMC draws detected")
    }

    #

    filt_vals <- obj$states()$filter_vals
    
    n_data <- dim(filt_vals)[1]
    n_states <- dim(filt_vals)[2]
    
    filt_vals_sorted <- apply(filt_vals,c(1,2),sort)
    filt_vals_sorted <- aperm(filt_vals_sorted,c(2,3,1))

    states_mean <- apply(filt_vals_sorted,c(1,2),mean)

    #

    upper_conf <- round(percentiles[3]*n_draws)
    mid_conf   <- round(percentiles[2]*n_draws)
    lower_conf <- round(percentiles[1]*n_draws)
    
    #

    plot_vals <- array(NA,dim=c(n_data,4,n_states))

    for (i in 1:n_states) {
        # Use the mean or selected percentile?
        if (use_mean == TRUE) {
            FDataTemp<-data.frame(filt_vals_sorted[,i,lower_conf],states_mean[,i],filt_vals_sorted[,i,upper_conf])
        } else {
            FDataTemp<-data.frame(filt_vals_sorted[,i,lower_conf],filt_vals_sorted[,i,mid_conf],filt_vals_sorted[,i,upper_conf])
        }

        FDataTemp <- as.matrix(FDataTemp)
        FDataTemp <- cbind(FDataTemp,1:n_data)

        #

        plot_vals[,,i] <- FDataTemp
    }

    #

    if (class(var_names) != "character") {
        var_names <- character(length=n_states)
        for (i in 1:n_states) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }
    
    #
    # plotting

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}

    #

    MR <- 0; MC <- 0
    plot_pages <- 1
    
    if(n_states < 4){
        MR <- n_states; MC <- 1
    }else if(n_states == 4){
        MR <- 2; MC <-2
    }else if(n_states > 4 && n_states < 7){
        MR <- 3; MC <- 2
    }else if(n_states > 6 && n_states < 10){
        MR <- 3; MC <- 3
    }else if(n_states > 9 && n_states < 13){
        MR <- 4; MC <- 3
    }else if(n_states > 12 && n_states < 25){
        MR <- 4; MC <- 3
        plot_pages <- 2
    }else if(n_states > 24 && n_states < 37){
        MR <- 4; MC <- 3
        plot_pages <- 3
    }else if(n_states > 36 && n_states < 49){
        MR <- 4; MC <- 3
        plot_pages <- 4
    }else if(n_states > 48 && n_states < 61){
        MR <- 4; MC <- 3
        plot_pages <- 5
    }else if(n_states > 60 && n_states < 73){
        MR <- 4; MC <- 3
        plot_pages <- 6
    }

    #

    state_count <- 1

    for (j in 1:plot_pages) {
        
        if(save==TRUE){
            if(class(dev.list()) != "NULL"){dev.off()}

            if(plot_pages==1){
                cairo_ps(filename="States.eps",height=height,width=width)
            }else{
                SaveState <- paste("States_",j,".eps",sep="")
                cairo_ps(filename=SaveState,height=height,width=width)
            }
        }
        
        #
        
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))
        
        #
        
        Time <- SL <- SU <- SM <- NULL # CRAN check workaround
        
        for (i in 1:MR) {
            for (k in 1:MC) {
                #
                if (state_count <= n_states) {
                    StateName <- var_names[state_count]
                    SDF <- plot_vals[,,state_count]
                    SDF <- data.frame(SDF)
                    colnames(SDF) <- c("SL","SM","SU","Time")

                    #

                    print(ggplot(data=SDF,aes(x=Time)) + xlab("Time") + ylab(paste(StateName)) + geom_ribbon(aes(ymin=SL,ymax=SU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_line(aes(y=SM),color="red",size=1) + theme(panel.background = element_rect(fill='white', colour='grey15')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                    
                    #
                    
                    state_count <- state_count + 1
                    Sys.sleep(0.6)
                    
                } else {state_count <- state_count + 1}
            }
        }
        
        if(save==TRUE){dev.off()}
    }
    #
    return=list(state_mean=states_mean,plot_vals=plot_vals)
}
