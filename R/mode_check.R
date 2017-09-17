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

mode_check.Rcpp_dsge_gensys <- function(obj,mode_vals=NULL,grid_size=201,scale_val=1,par_names=NULL,save=FALSE,height=13,width=13,...){
    .mode_check_int(obj,mode_vals,grid_size,scale_val,par_names,save,height,width) 
}

mode_check.Rcpp_dsge_uhlig <- function(obj,mode_vals=NULL,grid_size=201,scale_val=1,par_names=NULL,save=FALSE,height=13,width=13,...){
    .mode_check_int(obj,mode_vals,grid_size,scale_val,par_names,save,height,width) 
}

mode_check.Rcpp_dsgevar_gensys <- function(obj,mode_vals=NULL,grid_size=201,scale_val=1,par_names=NULL,save=FALSE,height=13,width=13,...){
    .mode_check_int(obj,mode_vals,grid_size,scale_val,par_names,save,height,width) 
}

mode_check.Rcpp_dsgevar_uhlig <- function(obj,mode_vals=NULL,grid_size=201,scale_val=1,par_names=NULL,save=FALSE,height=13,width=13,...){
    .mode_check_int(obj,mode_vals,grid_size,scale_val,par_names,save,height,width) 
}

#

.mode_check_int <- function(obj,mode_vals=NULL,grid_size=201,scale_val=1,par_names=NULL,save=FALSE,height=13,width=13){
    #
    
    n_param <- length(mode_vals)

    mode_check_vals <- obj$mode_check(mode_vals,grid_size,scale_val)$mode_check_vals

    #

    if (is.null(par_names)==TRUE) {
        par_names <- character(length=n_param)
        for(i in 1:n_param){  
            par_names[i] <- paste("parameter",i,sep="")
        }
    }

    #
    # Plot
    
    MR <- 0; MC <- 0
    plot_pages <- 1

        
    ParameterVals <- LogPosterior <- NULL # CRAN check workaround
    
    #
    
    if(n_param < 4){
        MR <- n_param; MC <- 1
    }else if(n_param == 4){
        MR <- 2; MC <-2
    }else if(n_param > 4 && n_param < 7){
        MR <- 3; MC <- 2
    }else if(n_param > 6 && n_param < 10){
        MR <- 3; MC <- 3
    }else if(n_param > 9 && n_param < 13){
        MR <- 4; MC <- 3
    }else if(n_param > 12 && n_param < 25){
        MR <- 4; MC <- 3
        plot_pages <- 2
    }else if(n_param > 24 && n_param < 37){
        MR <- 4; MC <- 3
        plot_pages <- 3
    }else if(n_param > 36 && n_param < 49){
        MR <- 4; MC <- 3
        plot_pages <- 4
    }else if(n_param > 48 && n_param < 61){
        MR <- 4; MC <- 3
        plot_pages <- 5
    }else if(n_param > 60 && n_param < 73){
        MR <- 4; MC <- 3
        plot_pages <- 6
    }
    
    #
    
    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    
    #
    
    param_count <- 1
    
    for (ii in 1:plot_pages) {
        if (save==TRUE) {
            if(class(dev.list()) != "NULL"){dev.off()}
            #
            SaveMode <- paste("DSGEModeCheck",as.character(ii),".eps",sep="")
            cairo_ps(filename=SaveMode,height=height,width=width)
        }
        
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))
        
        #
        
        for (i in 1:MR) {
            for (k in 1:MC) {
                #
                if (param_count <= (n_param)) {
                    ParamDF <- data.frame(mode_check_vals[,param_count,1],mode_check_vals[,param_count,2])
                    colnames(ParamDF) <- c("ParameterVals","LogPosterior")
                    
                    #
                    
                    print(ggplot(data=(ParamDF),aes(x=ParameterVals)) + xlab("") + ylab("Log Posterior") + geom_vline(xintercept=mode_vals[param_count],linetype = "longdash") + geom_line(aes(y=LogPosterior),color="springgreen4") + labs(title=par_names[param_count]) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                    
                    #
                    
                    param_count <- param_count + 1
                    Sys.sleep(0.3)
                    
                } else {param_count <- param_count + 1}
            }
        }
        
        if(save==TRUE){dev.off()}
    }
    #
}
