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

IRF.Rcpp_bvarm <- function(obj,periods=10,varnames=NULL,percentiles=c(.05,.50,.95),which_shock=NULL,which_response=NULL,save=TRUE,height=13,width=13,...)
{
    .irfvar(obj,periods,varnames,percentiles,which_shock,which_response,save,height,width)
}

IRF.Rcpp_bvars <- function(obj,periods=10,varnames=NULL,percentiles=c(.05,.50,.95),which_shock=NULL,which_response=NULL,save=TRUE,height=13,width=13,...)
{
    .irfvar(obj,periods,varnames,percentiles,which_shock,which_response,save,height,width)
}

IRF.Rcpp_bvarw <- function(obj,periods=10,varnames=NULL,percentiles=c(.05,.50,.95),which_shock=NULL,which_response=NULL,save=TRUE,height=13,width=13,...)
{
    .irfvar(obj,periods,varnames,percentiles,which_shock,which_response,save,height,width)
}

IRF.Rcpp_cvar <- function(obj,periods=10,varnames=NULL,percentiles=c(.05,.50,.95),which_shock=NULL,which_response=NULL,save=TRUE,height=13,width=13,...)
{
    .irfvar(obj,periods,varnames,percentiles,which_shock,which_response,save,height,width)
}

IRF.Rcpp_bvartvp <- function(obj,periods=10,which_irfs=NULL,varnames=NULL,percentiles=c(.05,.50,.95),which_shock=NULL,which_response=NULL,save=TRUE,height=13,width=13,...)
{
    .irfbvartvp(obj,periods,which_irfs,varnames,percentiles,which_shock,which_response,save,height,width)
}

IRF.Rcpp_gensys <- function(obj,periods=10,varnames=NULL,shocks_cov=NULL,save=FALSE,height=13,width=13,...){
    .irfdsge(obj,periods,varnames,shocks_cov,save,height,width)
}

IRF.Rcpp_uhlig <- function(obj,periods=10,varnames=NULL,shocks_cov=NULL,save=FALSE,height=13,width=13,...){
    .irfdsge(obj,periods,varnames,shocks_cov,save,height,width)
}

IRF.Rcpp_dsge_gensys <- function(obj,periods=10,obs_irfs=FALSE,varnames=NULL,percentiles=c(.05,.50,.95),save=TRUE,height=13,width=13,...){
    .irfedsge(obj,periods,obs_irfs,varnames,percentiles,save,height,width)
}

# IRF.EDSGE <- function(obj,observableIRFs=FALSE,varnames=NULL,percentiles=c(.05,.50,.95),save=TRUE,height=13,width=13,...){
#     .irfedsge(obj,observableIRFs,varnames,percentiles,save,height,width)
# }

IRF.Rcpp_dsgevar_gensys <- function(obj,periods,varnames=NULL,percentiles=c(.05,.50,.95),plot_comparison=TRUE,save=TRUE,height=13,width=13,...){
    .irfdsgevar(obj,periods,varnames,percentiles,plot_comparison,save,height,width)
}

# IRF.DSGEVAR <- function(obj,varnames=NULL,percentiles=c(.05,.50,.95),comparison=TRUE,save=TRUE,height=13,width=13,...){
#     .irfdsgevar(obj,varnames,percentiles,comparison,save,height,width)
# }


.irfvar <- function(obj, periods=10, varnames=NULL, percentiles=c(.05,.50,.95), which_shock=NULL, which_response=NULL, save=TRUE, height=13, width=13)
{
    #
    
    if (periods <= 0) {
        stop("error: need periods > 0")
    }

    #

    M <- obj$M
    n_draws <- dim(obj$beta_draws)[3]
    
    if (sum(dim(obj$irfs)) == 0) {
        obj$IRF(periods)
    }

    # put the IRFs in a tesseract-type format

    irf_temp <- obj$irfs # make a copy; much faster than accessing IRFs slice-by-slice in the loop below
    irf_tess <- array(NA,dim=c(M,M,periods,n_draws))

    for(i in 1:n_draws){
        # irf_tess[,,,i] <- obj$irfs[,,((i-1)*periods+1):(i*periods)]
        irf_tess[,,,i] <- irf_temp[,,((i-1)*periods+1):(i*periods)]
    }

    rm("irf_temp")

    irf_tess <- apply(irf_tess,c(3,1,2),sort)
    irf_tess <- aperm(irf_tess,c(2,3,1,4))

    irf_upper <- round(percentiles[3]*n_draws)
    irf_mid <- round(percentiles[2]*n_draws)
    irf_lower <- round(percentiles[1]*n_draws)

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}

    if(class(dev.list()) != "NULL"){dev.off()}
    
    #

    if (is.null(which_shock)) {
        which_shock <- 1:M
    }

    if (is.null(which_response)) {
        which_response <- 1:M
    }

    n_response <- length(which_response)
    n_shocks   <- length(which_shock)

    if (class(varnames) != "character") {
        varnames <- character(length=M)
        for (i in 1:M) {  
            varnames[i] <- paste("VAR",i,sep="")
        }
    }

    #

    irf_plot <- array(NA,dim=c(periods,4,M,M))

    for(i in 1:M){
        for(k in 1:M){
            IRFPData <- data.frame(irf_tess[,k,irf_lower,i],irf_tess[,k,irf_mid,i],irf_tess[,k,irf_upper,i],1:(periods))
            IRFPData <- as.matrix(IRFPData)
            irf_plot[,,k,i] <- IRFPData
        }
    }

    #
    # plot IRFs

    IRFL <- IRFM <- IRFU <- Time <- NULL # CRAN check workaround
    
    if (n_response == M && n_shocks == M) {

        if(class(dev.list()) != "NULL"){dev.off()}
        if(save==TRUE){cairo_ps(filename="IRFs.eps",height=height,width=width)}
        
        pushViewport(viewport(layout=grid.layout(M,M)))
        
        for(i in 1:M){
            for(k in 1:M){
                NameResponse <- varnames[k]
                NameImpulse  <- varnames[i]
                
                IRFDF <- irf_plot[,,k,i]
                IRFDF <- data.frame(IRFDF)
                colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")
                #
                gg1 <- ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste("Shock from ",NameImpulse," to", NameResponse)) 
                gg2 <- gg1 + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="red",size=2) 
                gg3 <- gg2 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
                print(gg3,vp = vplayout(i,k))
                #
                Sys.sleep(0.3)
            }
        }
        if (save==TRUE) {dev.off()}
    } else {
        
        if (n_response < 4) {
            MR <- n_response; MC <- 1
        } else if (n_response == 4) {
            MR <- 2; MC <-2
        } else if (n_response > 4 && n_response < 7) {
            MR <- 3; MC <- 2
        } else if (n_response > 6 && n_response < 10) {
            MR <- 3; MC <- 3
        } else if (n_response > 9 && n_response < 13) {
            MR <- 4; MC <- 3
        }else if (n_response > 12 && n_response < 17) {
            MR <- 4; MC <- 4
        }else if (n_response > 17 && n_response < 21) {
            MR <- 5; MC <- 4
        }else if (n_response > 20 && n_response < 26) {
            MR <- 5; MC <- 5
        } else if (n_response > 25 && n_response < 31) {
            MR <- 5; MC <- 6
        } else if (n_response > 30 && n_response < 37) {
            MR <- 6; MC <- 6
        } else {
            stop("You have too many IRFs to plot!")
        }

        #
        
        for (i in which_shock) {
            plot_ind_r <- 1
            plot_ind_c <- 1
            
            if (save==TRUE) {
                if (n_shocks==1) {
                    cairo_ps(filename="IRFs.eps",height=height,width=width)
                } else {
                    SaveIRF <- paste(varnames[i],"_Shock",".eps",sep="")
                    cairo_ps(filename=SaveIRF,height=height,width=width)
                }
            }

            grid.newpage()
            pushViewport(viewport(layout=grid.layout(MR,MC)))
            
            for(k in which_response){
                NameResponse <- varnames[k]
                NameImpulse  <- varnames[i]
                
                IRFDF <- irf_plot[,,k,i]
                IRFDF <- data.frame(IRFDF)

                colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")

                #

                gg1 <- ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste("Shock from ", NameImpulse," to", NameResponse)) 
                gg2 <- gg1 + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="red",size=2) 
                gg3 <- gg2 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))

                print(gg3,vp = vplayout(plot_ind_r,plot_ind_c))
                
                Sys.sleep(0.3)

                #

                plot_ind_c <- plot_ind_c + 1

                if (plot_ind_c > MC) {
                    plot_ind_c <- 1
                    plot_ind_r <- plot_ind_r + 1
                }
            }
        }

        if(save==TRUE){dev.off()}
    }
    #
    return=list(IRFs=irf_plot)
}

.irfbvartvp <- function(obj, periods=10, which_irfs=NULL, varnames=NULL, percentiles=c(.05,.50,.95), which_shock=NULL, which_response=NULL, save=TRUE, height=13, width=13)
{
    #
    
    if (periods <= 0) {
        stop("error: need periods > 0")
    }

    #

    M <- obj$M
    n_draws <- dim(obj$alpha_draws)[3]

    irf_upper <- round(percentiles[3]*n_draws)
    irf_mid <- round(percentiles[2]*n_draws)
    irf_lower <- round(percentiles[1]*n_draws)

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}

    if(class(dev.list()) != "NULL"){dev.off()}
    
    #

    if (is.null(which_shock)) {
        which_shock <- 1:M
    }

    if (is.null(which_response)) {
        which_response <- 1:M
    }

    n_response <- length(which_response)
    n_shocks   <- length(which_shock)

    if (class(varnames) != "character") {
        varnames <- character(length=M)
        for (i in 1:M) {  
            varnames[i] <- paste("VAR",i,sep="")
        }
    }

    #
    # loop

    if (is.null(which_irfs)) {
        which_irfs <- dim(obj$beta_draws)[2]
    }

    n_plots <- length(which_irfs)

    for (pp in 1:n_plots) {

        obj$IRF(periods,which_irfs[pp]-1) # minus 1 to account for zero-index

        # put the IRFs in a tesseract-type format

        irf_temp <- obj$irfs
        irf_tess <- array(NA,dim=c(M,M,periods,n_draws))

        for(i in 1:n_draws){
            irf_tess[,,,i] <- irf_temp[,,((i-1)*periods+1):(i*periods)]
        }

        rm("irf_temp")

        irf_tess <- apply(irf_tess,c(3,1,2),sort)
        irf_tess <- aperm(irf_tess,c(2,3,1,4))
        
        irf_plot <- array(NA,dim=c(periods,4,M,M))

        for(i in 1:M){
            for(k in 1:M){
                IRFPData <- data.frame(irf_tess[,k,irf_lower,i],irf_tess[,k,irf_mid,i],irf_tess[,k,irf_upper,i],1:(periods))
                IRFPData <- as.matrix(IRFPData)
                irf_plot[,,k,i] <- IRFPData
            }
        }

        #
        # plot IRFs

        IRFL <- IRFM <- IRFU <- Time <- NULL # CRAN check workaround
        
        if (n_response == M && n_shocks == M) {

            if(save==TRUE){cairo_ps(filename="IRFs.eps",height=height,width=width)}
            
            grid.newpage()
            pushViewport(viewport(layout=grid.layout(M,M)))
            
            for(i in 1:M){
                for(k in 1:M){
                    NameResponse <- varnames[k]
                    NameImpulse  <- varnames[i]
                    
                    IRFDF <- irf_plot[,,k,i]
                    IRFDF <- data.frame(IRFDF)
                    colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")
                    #
                    gg1 <- ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste("Shock from ",NameImpulse," to", NameResponse)) 
                    gg2 <- gg1 + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="red",size=2) 
                    gg3 <- gg2 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
                    print(gg3,vp = vplayout(i,k))
                    #
                    Sys.sleep(0.3)
                }
            }
            if (save==TRUE) {dev.off()}
        } else {
            
            if (n_response < 4) {
                MR <- n_response; MC <- 1
            } else if (n_response == 4) {
                MR <- 2; MC <-2
            } else if (n_response > 4 && n_response < 7) {
                MR <- 3; MC <- 2
            } else if (n_response > 6 && n_response < 10) {
                MR <- 3; MC <- 3
            } else if (n_response > 9 && n_response < 13) {
                MR <- 4; MC <- 3
            }else if (n_response > 12 && n_response < 17) {
                MR <- 4; MC <- 4
            }else if (n_response > 17 && n_response < 21) {
                MR <- 5; MC <- 4
            }else if (n_response > 20 && n_response < 26) {
                MR <- 5; MC <- 5
            } else if (n_response > 25 && n_response < 31) {
                MR <- 5; MC <- 6
            } else if (n_response > 30 && n_response < 37) {
                MR <- 6; MC <- 6
            } else {
                stop("You have too many IRFs to plot!")
            }

            #
            
            for (i in which_shock) {
                plot_ind_r <- 1
                plot_ind_c <- 1
                
                if (save==TRUE) {
                    if (n_shocks==1) {
                        cairo_ps(filename="IRFs.eps",height=height,width=width)
                    } else {
                        SaveIRF <- paste(varnames[i],"_Shock",".eps",sep="")
                        cairo_ps(filename=SaveIRF,height=height,width=width)
                    }
                }

                grid.newpage()
                pushViewport(viewport(layout=grid.layout(MR,MC)))
                
                for(k in which_response){
                    NameResponse <- varnames[k]
                    NameImpulse  <- varnames[i]
                    
                    IRFDF <- irf_plot[,,k,i]
                    IRFDF <- data.frame(IRFDF)

                    colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")

                    #

                    gg1 <- ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste("Shock from ", NameImpulse," to", NameResponse)) 
                    gg2 <- gg1 + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="red",size=2) 
                    gg3 <- gg2 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))

                    print(gg3,vp = vplayout(plot_ind_r,plot_ind_c))
                    
                    Sys.sleep(0.3)

                    #

                    plot_ind_c <- plot_ind_c + 1

                    if (plot_ind_c > MC) {
                        plot_ind_c <- 1
                        plot_ind_r <- plot_ind_r + 1
                    }
                }
            }

            if(save==TRUE){dev.off()}
        }
    }
    #
    return=list(IRFs=irf_plot)
}

#
# DSGE

.irfdsge <- function(obj,periods=10,varnames=NULL,shocks_cov=NULL,save=FALSE,height=13,width=13)
{
    
    if (periods <= 0) {
        stop("error: need periods > 0")
    }

    #

    irfs <- obj$IRF(periods)$irf_vals

    irfs <- round(irfs,10)

    M <- dim(irfs)[2]
    n_shocks <- dim(irfs)[3]

    if (class(varnames) != "character") {
        varnames <- character(length=M)
        for (i in 1:M) {  
            varnames[i] <- paste("VAR",i,sep="")
        }
    }
    
    #
    
    irfs_plot <- irfs
    drop_ind <- (M - n_shocks + 1):M
    
    if (n_shocks > 1) { # drop other shocks
        irfs_plot <- array(0, dim=c(periods, M - n_shocks + 1, n_shocks))
        
        for (j in 1:n_shocks) {
            irfs_plot[,,j] <- irfs[,-drop_ind[-(j)],j]
        }
    }
    
    #
    # Plot
    
    IRFM <- Time <- NULL # CRAN check workaround
    
    n_response <- M - n_shocks + 1
    MR <- 0; MC <- 0
    
    if (n_response < 4) {
        MR <- n_response; MC <- 1
    } else if (n_response == 4) {
        MR <- 2; MC <-2
    } else if (n_response > 4 && n_response < 7) {
        MR <- 3; MC <- 2
    } else if (n_response > 6 && n_response < 10) {
        MR <- 3; MC <- 3
    } else if (n_response > 9 && n_response < 13) {
        MR <- 4; MC <- 3
    } else if (n_response > 12 && n_response < 17) {
        MR <- 4; MC <- 4
    } else if (n_response > 17 && n_response < 21) {
        MR <- 5; MC <- 4
    } else if (n_response > 20 && n_response < 26) {
        MR <- 5; MC <- 5
    } else if (n_response > 25 && n_response < 31) {
        MR <- 5; MC <- 6
    } else if (n_response > 30 && n_response < 37) {
        MR <- 6; MC <- 6
    } else {
        stop("You have too many IRFs to plot!")
    }
    
    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}

    if(class(dev.list()) != "NULL"){dev.off()}
    
    #

    varnames2 <- varnames

    for (j in 1:n_shocks) {
        
        if (save==TRUE) {
            if (n_shocks==1) {
                cairo_ps(filename="DSGEIRFs.eps",height=height,width=width)
            } else {
                SaveIRF <- paste(varnames[M-n_shocks+j],"_Shock",".eps",sep="")
                cairo_ps(filename=SaveIRF,height=height,width=width)
            }
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))

        if (n_shocks > 1) {
            varnames2 <- varnames[-drop_ind[-(j)]]
        } else {
            varnames2 <- varnames
        }

        #

        irf_plot_count <- 1

        for (i in 1:MR) {
            for (k in 1:MC) {
                
                if (irf_plot_count <= (M - n_shocks + 1)) {
                    IRFDF <- data.frame(irfs_plot[,irf_plot_count,j],1:periods)
                    colnames(IRFDF) <- c("IRFM","Time")
                    
                    print(ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(varnames2[irf_plot_count])) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="darkslateblue",size=2) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                    
                    #

                    irf_plot_count <- irf_plot_count + 1
                    
                    Sys.sleep(0.3)
                } else {
                    irf_plot_count <- irf_plot_count + 1
                }
            }
        }

        if(save==TRUE){dev.off()}
    }
    #
    return=list(irfs=irfs)
}

.irfedsge <- function(obj,periods=10,obs_irfs=FALSE,varnames=NULL,percentiles=c(.05,.50,.95),save=TRUE,height=13,width=13)
{    
    if (periods <= 0) {
        stop("error: need periods > 0")
    }

    #

    n_draws <- dim(obj$dsge_draws)[1]

    if (n_draws <= 0) {
        stop("error: no MCMC draws detected")
    }

    irfs <- obj$IRF(periods)$irf_vals # deep copy is needed
    irfs <- round(irfs,10)

    M <- dim(irfs)[2]
    n_shocks <- dim(irfs)[3] / n_draws

    n_response <- M - n_shocks + 1

    #

    irfs_i <- array(0, dim=c(periods, M, n_shocks))
    irfs_tess <- array(0, dim=c(periods, n_response, n_shocks, n_draws))

    drop_ind <- n_response:M
    
    for (i in 1:n_draws) {
        irfs_i <- irfs[,,((i-1)*n_shocks + 1):(i*n_shocks)]

        if (n_shocks > 1) {
            for (j in 1:n_shocks) {
                irfs_tess[,,j,i] <- irfs_i[,-drop_ind[-(j)],j]
            }
        } else {
            irfs_tess[,,,i] <- irfs_i
        }
    }

    irfs_tess <- aperm(irfs_tess,c(2,3,1,4))
    irfs_tess <- apply(irfs_tess,c(3,1,2),sort)

    #

    irf_upper <- round(percentiles[3]*n_draws)
    irf_mid <- round(percentiles[2]*n_draws)
    irf_lower <- round(percentiles[1]*n_draws)

    #

    irfs_plot <- array(NA,dim=c(periods,4,n_response,n_shocks))

    for (i in 1:n_shocks) {
        for (k in 1:n_response) {
            irfs_plot_temp <- data.frame(c(irfs_tess[irf_lower,,k,i]),c(irfs_tess[irf_mid,,k,i]),c(irfs_tess[irf_upper,,k,i]),1:(periods))
            irfs_plot[,,k,i] <- data.matrix(irfs_plot_temp)
        }
    }

    #

    if (class(varnames) != "character") {
        varnames <- character(length=M)
        for (i in 1:M) {  
            varnames[i] <- paste("VAR",i,sep="")
        }
    }

    #

    IRFM <- IRFU <- IRFL <- Time <- NULL # CRAN check workaround

    MR <- 0; MC <- 0
    if(n_response < 4){
        MR <- n_response; MC <- 1
    }else if(n_response == 4){
        MR <- 2; MC <-2
    }else if(n_response > 4 && n_response < 7){
        MR <- 3; MC <- 2
    }else if(n_response > 6 && n_response < 10){
        MR <- 3; MC <- 3
    }else if(n_response > 9 && n_response < 13){
        MR <- 4; MC <- 3
    }else if(n_response > 12 && n_response < 17){
        MR <- 4; MC <- 4
    }else if(n_response > 17 && n_response < 21){
        MR <- 5; MC <- 4
    }else if(n_response > 20 && n_response < 26){
        MR <- 5; MC <- 5
    }else if(n_response > 25 && n_response < 31){
        MR <- 5; MC <- 6
    }else if(n_response > 30 && n_response < 37){
        MR <- 6; MC <- 6
    }else{
        stop("You have too many IRFs to plot!")
    }

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}

    if(class(dev.list()) != "NULL"){dev.off()}
    
    #

    varnames2 <- varnames

    for (j in 1:n_shocks) {
        #
        if (save==TRUE) {
            if (n_shocks==1) {
                cairo_ps(filename="DSGEBIRFs.eps",height=height,width=width)
            } else {
                SaveIRF <- paste(varnames[M-n_shocks+j],"_Shock",".eps",sep="")
                #
                cairo_ps(filename=SaveIRF,height=height,width=width)
            }
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))

        #

        if (n_shocks > 1) {
            varnames2 <- varnames[-drop_ind[-(j)]]
        } else {
            varnames2 <- varnames
        }

        #

        irf_plot_count <- 1

        for (i in 1:MR) {
            for (k in 1:MC) {
                #
                if (irf_plot_count <= n_response) {
                    IRFDF <- data.frame(irfs_plot[,,irf_plot_count,j])
                    colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")

                    #

                    gg1 <- ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(varnames2[irf_plot_count])) 
                    gg2 <- gg1 + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="darkslateblue",size=2) 
                    gg3 <- gg2 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
                    print(gg3,vp = vplayout(i,k))
                    
                    #

                    irf_plot_count <- irf_plot_count + 1
                    
                    Sys.sleep(0.3)
                } else {
                    irf_plot_count <- irf_plot_count + 1
                }
            }
        }

        if(save==TRUE){dev.off()}
    }
    #
    return=list()
}

.irfdsgevar <- function(obj,periods,varnames=NULL,percentiles=c(.05,.50,.95),comparison_plot=TRUE,save=TRUE,height=13,width=13)
{   
    if (periods <= 0) {
        stop("error: need periods > 0")
    }

    #

    n_draws <- dim(obj$beta_draws)[3]

    if (n_draws <= 0) {
        stop("error: no MCMC draws detected")
    }

    irfs <- obj$IRF(periods)$irf_vals # deep copy is needed
    irfs <- round(irfs,10)

    M <- dim(irfs)[2] 
    n_shocks <- M

    n_response <- M 

    #

    irfs_tess <- array(0, dim=c(M, M, periods, n_draws))
    
    for (i in 1:n_draws) {
        irfs_tess[,,,i] <- irfs[,,((i-1)*periods + 1):(i*periods)]
    }

    irfs_tess <- apply(irfs_tess,c(3,1,2),sort)
    #irfs_tess <- aperm(irfs_tess,c(2,3,1,4))
    #irfs_tess <- aperm(irfs_tess,c(1,2,4,3))
    #irfs_tess <- aperm(irfs_tess,c(4,1,2,3))

    #

    irf_upper <- round(percentiles[3]*n_draws)
    irf_mid <- round(percentiles[2]*n_draws)
    irf_lower <- round(percentiles[1]*n_draws)

    #

    irfs_plot <- array(NA,dim=c(periods,4,n_response,n_shocks))

    for (i in 1:M) {
        for (k in 1:M) {
            irfs_plot_temp <- data.frame(c(irfs_tess[irf_lower,,k,i]),c(irfs_tess[irf_mid,,k,i]),c(irfs_tess[irf_upper,,k,i]),1:(periods))
            irfs_plot[,,k,i] <- data.matrix(irfs_plot_temp)
        }
    }
    
    #

    if (class(varnames) != "character") {
        varnames <- character(length=M)
        for (i in 1:M) {  
            varnames[i] <- paste("VAR",i,sep="")
        }
    }

    #

    IRFM <- IRFU <- IRFL <- Time <- NULL # CRAN check workaround

    MR <- 0; MC <- 0
    if(n_response < 4){
        MR <- n_response; MC <- 1
    }else if(n_response == 4){
        MR <- 2; MC <-2
    }else if(n_response > 4 && n_response < 7){
        MR <- 3; MC <- 2
    }else if(n_response > 6 && n_response < 10){
        MR <- 3; MC <- 3
    }else if(n_response > 9 && n_response < 13){
        MR <- 4; MC <- 3
    }else if(n_response > 12 && n_response < 17){
        MR <- 4; MC <- 4
    }else if(n_response > 17 && n_response < 21){
        MR <- 5; MC <- 4
    }else if(n_response > 20 && n_response < 26){
        MR <- 5; MC <- 5
    }else if(n_response > 25 && n_response < 31){
        MR <- 5; MC <- 6
    }else if(n_response > 30 && n_response < 37){
        MR <- 6; MC <- 6
    }else{
        stop("You have too many IRFs to plot!")
    }

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}

    if(class(dev.list()) != "NULL"){dev.off()}
    
    #

    for (j in 1:n_shocks) {
        #
        if (save==TRUE) {
            if (n_shocks==1) {
                cairo_ps(filename="DSGEVARIRFs.eps",height=height,width=width)
            } else {
                SaveIRF <- paste("Shock_",j,".eps",sep="")
                cairo_ps(filename=SaveIRF,height=height,width=width)
            }
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))

        #

        irf_plot_count <- 1

        for (i in 1:MR) {
            for (k in 1:MC) {
                #
                if (irf_plot_count <= n_response) {
                    
                    IRFDF <- data.frame(irfs_plot[,,irf_plot_count,j])
                    colnames(IRFDF) <- c("VARIRFL","VARIRFM","VARIRFU","Time")

                    #
                    
                    print(ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(varnames[irf_plot_count])) + geom_hline(yintercept=0) + geom_line(aes(y=VARIRFM),color="darkgreen",size=2) + geom_line(aes(y=VARIRFL),color="darkgreen",size=1,linetype=4) + geom_line(aes(y=VARIRFU),color="darkgreen",size=1,linetype=4) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))

                    #
                    
                    irf_plot_count <- irf_plot_count + 1
                    Sys.sleep(0.3)
                    
                } else {
                    irf_plot_count <- irf_plot_count + 1
                }
            }
        }

        if(save==TRUE){dev.off()}
    }
    #
    return=list()
}

# old

# .irfbvartvp <- function(obj, which_irfs=NULL, varnames=NULL, percentiles=c(.05,.50,.95), save=FALSE, height=13, width=13){
#     IRFs <- obj$IRFs
#     irf.periods <- as.numeric(dim(IRFs)[1])
#     M <- as.numeric(dim(IRFs)[3])
#     irf.points <- obj$irf.points
#     mydata <- obj$data
#     tau <- obj$tau
#     #
#     if(class(dev.list()) != "NULL"){dev.off()}
#     #
#     vplayout<-function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
#     #
#     total <- as.numeric(dim(IRFs)[5])
#     IRFUpper <- round(percentiles[3]*total)
#     IRFMid <- round(percentiles[2]*total)
#     IRFLower <- round(percentiles[1]*total)
#     #
#     if(class(whichirfs)!="NULL"){
#         imp.pl <- as.numeric(length(whichirfs))
#         whichirfs <- whichirfs - tau
#     }else{
#         imp.pl <- 1:as.numeric(dim(IRFs)[4])
#         whichirfs <- imp.pl
#     }
#     #
#     IRFPlot<-array(NA,dim=c(irf.periods,4,M,M,as.numeric(dim(IRFs)[4])))
#     for(ik in imp.pl){
#         iki <- whichirfs[ik]
#         for(i in 1:M){ #shock number
#             for(k in 1:M){ #variable number
#                 IRFPData<-data.frame(IRFs[,k,i,iki,IRFLower],IRFs[,k,i,iki,IRFMid],IRFs[,k,i,iki,IRFUpper],1:(irf.periods))
#                 IRFPData<-as.matrix(IRFPData)
#                 IRFPlot[,,k,i,ik]<-IRFPData
#             }
#         }
#     }
#     #
#     IRFL <- IRFM <- IRFU <- Time <- NULL # CRAN check workaround
#     value <- variable <- NULL
#     #
#     for(ik in 1:as.numeric(dim(IRFs)[4])){
#         SaveName <- paste(as.character(irf.points[ik]),".eps",sep="")
#         if(save==TRUE){cairo_ps(filename=SaveName,height=height,width=width)}
#         #
#         grid.newpage()
#         pushViewport(viewport(layout=grid.layout(M,M)))
#         #
#         for(i in 1:M){
#             for(k in 1:M){
#                 NameResponse <- colnames(mydata)[k]
#                 NameImpulse <- colnames(mydata)[i]
#                 #
#                 IRFDF <- IRFPlot[,,k,i,ik]
#                 IRFDF <- data.frame(IRFDF)
#                 colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")
#                 #
#                 print(ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste("Shock from ",NameImpulse," to", NameResponse)) + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="red",size=2),vp = vplayout(i,k))
#                 #
#                 Sys.sleep(0.3)
#             }
#         }
#         if(save==TRUE){dev.off()}
#     }
#     #
#     SaveCompName <- "IRFComparison.eps"
#     if(save==TRUE){cairo_ps(filename=SaveCompName,height=height,width=width)}
#     grid.newpage()
#     pushViewport(viewport(layout=grid.layout(M,M)))
#     for(i in 1:M){
#         for(k in 1:M){
#             NameResponse <- colnames(mydata)[k]
#             NameImpulse <- colnames(mydata)[i]
#             #
#             IRFDF <-data.frame(IRFPlot[,2,k,i,],1:irf.periods)
#             nnames <- c(as.character(irf.points),"irf.periods")
#             colnames(IRFDF) <- nnames
#             #
#             colIRFDF<-as.numeric(ncol(IRFDF))
#             rowIRFDF<-as.numeric(nrow(IRFDF))
#             newIRFDF <- matrix(NA,nrow=(nrow(IRFDF)*(ncol(IRFDF)-1)),ncol=3)
#             newIRFDF <- data.frame(newIRFDF)
#             for(j in 1:(ncol(IRFDF)-1)){
#                 newIRFDF[(((j-1)*rowIRFDF)+1):(j*rowIRFDF),1] <- 1:irf.periods
#                 newIRFDF[(((j-1)*rowIRFDF)+1):(j*rowIRFDF),2] <- rep(colnames(IRFDF)[j],rowIRFDF)
#             }
#             #
#             newIRFDF[,3] <- c(IRFPlot[,2,k,i,])
#             #
#             IRFDF <- newIRFDF
#             colnames(IRFDF) <- c("irf.periods","variable","value")
#             #
#             if(i == 1 & k == 1){
#                 print(ggplot(data=IRFDF,aes(x=irf.periods,y=value,colour=variable)) + geom_line() + xlab("") + ylab(paste("Shock from ",NameImpulse," to", NameResponse)) + geom_hline(yintercept=0) + scale_colour_hue("") + theme(legend.position="top",legend.direction="horizontal", legend.background=element_blank()),vp = vplayout(i,k))
#             }else{
#                 print(ggplot(data=IRFDF,aes(x=irf.periods,y=value,colour=variable)) + geom_line() + xlab("") + ylab(paste("Shock from ",NameImpulse," to", NameResponse)) + geom_hline(yintercept=0) + theme(legend.position="none"),vp = vplayout(i,k))
#             }
#             #
#             Sys.sleep(0.3)
#         }
#     }
#     if(save==TRUE){dev.off()}
#     #
# }
