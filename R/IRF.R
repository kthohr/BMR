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

# IRF.gensys <- function(obj,shocks,irf.periods=20,varnames=NULL,plot=TRUE,save=FALSE,height=13,width=13,...){
#     .irfsdsge(obj,shocks,irf.periods,varnames,plot,save,height,width)
# }

# IRF.uhlig <- function(obj,shocks,irf.periods=20,varnames=NULL,plot=TRUE,save=FALSE,height=13,width=13,...){
#     .irfsdsge(obj,shocks,irf.periods,varnames,plot,save,height,width)
# }

# IRF.EDSGE <- function(obj,observableIRFs=FALSE,varnames=NULL,percentiles=c(.05,.50,.95),save=TRUE,height=13,width=13,...){
#     .irfedsge(obj,observableIRFs,varnames,percentiles,save,height,width)
# }

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

.irfdsge <- function(obj,periods=10,varnames=NULL,shocks_cov=NULL,save=FALSE,height=13,width=13)
{
    
    if (periods <= 0) {
        stop("error: need periods > 0")
    }
    
    if (!is.null(shocks_cov)) {
        obj$shocks_cov <- shocks_cov
    }

    irfs <- obj$IRF(periods)$irf_vals

    irfs[abs(irfs) < 1e-14] <- 0

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

.irfedsge <- function(obj,observableIRFs=FALSE,varnames=NULL,percentiles=c(.05,.50,.95),save=TRUE,height=13,width=13){
    #
    if(nrow(obj$Parameters)==0){
        stop("no MCMC draws detected.\n",call.=FALSE)
    }
    #
    IRFs <- round(obj$IRFs,10)
    irf.periods <- as.numeric(dim(IRFs)[1])
    nResp <- as.numeric(dim(IRFs)[2])
    n_shocks <- as.numeric(dim(IRFs)[3])
    keep <- dim(IRFs)[4]
    #
    ObserveMat <- obj$ObserveMat
    #
    IRFObservable <- array(0,dim=c(irf.periods,ncol(ObserveMat),n_shocks,keep))
    if(observableIRFs == TRUE){
        if(ncol(ObserveMat)==1){
            for(i in 1:keep){
                for(j in 1:n_shocks){
                    IRFObservable[,1,j,i] <- IRFs[,,j,i]%*%ObserveMat
                }
            }
        }else{
            for(i in 1:keep){
                for(j in 1:n_shocks){
                    IRFObservable[,,j,i] <- IRFs[,,j,i]%*%ObserveMat
                }
            }
        }
        IRFs <- IRFObservable
        nResp <- as.numeric(ncol(ObserveMat))
    }
    #
    ImpSorted<-aperm(IRFs,c(2,3,1,4))
    #
    IRFs<-apply(ImpSorted,c(3,1,2),sort)
    #
    IRFUpper <- round(percentiles[3]*keep)
    IRFMid <- round(percentiles[2]*keep)
    IRFLower <- round(percentiles[1]*keep)
    #
    IRFPlot <- array(NA,dim=c(irf.periods,4,nResp,n_shocks))
    for(i in 1:n_shocks){
        for(k in 1:nResp){
            IRFPData <- data.frame(c(IRFs[IRFLower,,k,i]),c(IRFs[IRFMid,,k,i]),c(IRFs[IRFUpper,,k,i]),1:(irf.periods))
            IRFPData <- as.matrix(IRFPData)
            IRFPlot[,,k,i] <- IRFPData
        }
    }
    #
    IRFsRet <- IRFPlot
    #
    if(class(varnames) != "character"){
        varnames <- character(length=nResp)
        for(i in 1:nResp){  
            varnames[i] <- paste("VAR",i,sep="")
        }
    }
    #
    if(n_shocks > 1 && observableIRFs!=TRUE){
        IRFPlot <- array(0,dim=c(irf.periods,4,nResp-n_shocks+1,n_shocks))
        DropCount <- (nResp-n_shocks+1):nResp
        for(j in 1:n_shocks){
            IRFPlot[,,,j] <- IRFsRet[,,-DropCount[-(j)],j]
        }
    }
    #
    n_response <- nResp - n_shocks + 1
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
    if(class(dev.list()) != "NULL"){dev.off()}
    #
    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    #
    for(j in 1:n_shocks){
        #
        if(save==TRUE){
            if(n_shocks==1){
                cairo_ps(filename="DSGEBIRFs.eps",height=height,width=width)
            }else{
                SaveIRF <- paste(varnames[nResp-n_shocks+j],"_Shock",".eps",sep="")
                #
                cairo_ps(filename=SaveIRF,height=height,width=width)
            }
        }
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))
        #
        IRFM <- IRFU <- IRFL <- Time <- NULL # CRAN check workaround
        #
        if(n_shocks > 1 && observableIRFs!=TRUE){
            varnames2 <- varnames[-DropCount[-(j)]]
        }else{varnames2 <- varnames}
        #
        irf_plot_count <- 1
        IRFO <- 1
        if(observableIRFs==TRUE){IRFO <- nResp}
        #
        for(i in 1:MR){
            for(k in 1:MC){
                #
                if(irf_plot_count <= (nResp - n_shocks + IRFO)){
                    IRFDF <- data.frame(IRFPlot[,,irf_plot_count,j])
                    colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")
                    #
                    gg1 <- ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(varnames2[irf_plot_count])) 
                    gg2 <- gg1 + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="darkslateblue",size=2) 
                    gg3 <- gg2 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
                    print(gg3,vp = vplayout(i,k))
                    #
                    irf_plot_count <- irf_plot_count + 1
                    #
                    Sys.sleep(0.3)
                }else{irf_plot_count <- irf_plot_count + 1}
            }
        }
        if(save==TRUE){dev.off()}
    }
    #
    return=list(IRFs=IRFsRet)
}

.irfdsgevar <- function(obj,varnames=NULL,percentiles=c(.05,.50,.95),comparison=TRUE,save=TRUE,height=13,width=13){
    #
    if(nrow(obj$Parameters)==0){
        stop("no MCMC draws detected.\n",call.=FALSE)
    }
    #
    observableIRFs <- TRUE
    #
    IRFDs <- round(obj$DSGEIRFs,10)
    irf.periods <- as.numeric(dim(IRFDs)[1])
    nResp <- as.numeric(dim(IRFDs)[2])
    n_shocks <- as.numeric(dim(IRFDs)[3])
    keep <- dim(IRFDs)[4]
    #
    ObserveMat <- obj$ObserveMat
    #
    IRFObservable <- array(0,dim=c(irf.periods,ncol(ObserveMat),n_shocks,keep))
    if(observableIRFs == TRUE){
        if(ncol(ObserveMat)==1){
            for(i in 1:keep){
                for(j in 1:n_shocks){
                    IRFObservable[,1,j,i] <- IRFDs[,,j,i]%*%ObserveMat
                }
            }
        }else{
            for(i in 1:keep){
                for(j in 1:n_shocks){
                    IRFObservable[,,j,i] <- IRFDs[,,j,i]%*%ObserveMat
                }
            }
        }
        IRFDs <- IRFObservable
        nResp <- as.numeric(ncol(ObserveMat))
    }
    #
    ImpSorted<-aperm(IRFDs,c(2,3,1,4))
    #
    IRFDs<-apply(ImpSorted,c(3,1,2),sort)
    #
    IRFDVs <-aperm(obj$DSGEVARIRFs,c(4,1,2,3))
    #
    IRFUpper <- round(percentiles[3]*keep)
    IRFMid <- round(percentiles[2]*keep)
    IRFLower <- round(percentiles[1]*keep)
    #
    IRFDPlot <- array(NA,dim=c(irf.periods,4,nResp,n_shocks))
    IRFDVPlot <- array(NA,dim=c(irf.periods,4,nResp,n_shocks))
    for(i in 1:n_shocks){
        for(k in 1:nResp){
            IRFPData <- data.frame(c(IRFDs[IRFLower,,k,i]),c(IRFDs[IRFMid,,k,i]),c(IRFDs[IRFUpper,,k,i]),1:(irf.periods))
            IRFPData <- as.matrix(IRFPData)
            IRFDPlot[,,k,i] <- IRFPData
            #
            IRFPData <- data.frame(c(IRFDVs[IRFLower,,k,i]),c(IRFDVs[IRFMid,,k,i]),c(IRFDVs[IRFUpper,,k,i]),1:(irf.periods))
            IRFPData <- as.matrix(IRFPData)
            IRFDVPlot[,,k,i] <- IRFPData
        }
    }
    #
    IRFsRet <- IRFDPlot
    #
    if(class(varnames) != "character"){
        varnames <- character(length=nResp)
        for(i in 1:nResp){  
            varnames[i] <- paste("VAR",i,sep="")
        }
    }
    #
    if(n_shocks > 1 && observableIRFs!=TRUE){
        IRFDPlot <- array(0,dim=c(irf.periods,4,nResp-n_shocks+1,n_shocks))
        DropCount <- (nResp-n_shocks+1):nResp
        for(j in 1:n_shocks){
            IRFDPlot[,,,j] <- IRFsRet[,,-DropCount[-(j)],j]
        }
    }
    #
    MR <- 0; MC <- 0
    if(nResp < 4){
        MR <- nResp; MC <- 1
    }else if(nResp == 4){
        MR <- 2; MC <-2
    }else if(nResp > 4 && nResp < 7){
        MR <- 3; MC <- 2
    }else if(nResp > 6 && nResp < 10){
        MR <- 3; MC <- 3
    }else if(nResp > 9 && nResp < 13){
        MR <- 4; MC <- 3
    }else if(nResp > 12 && nResp < 17){
        MR <- 4; MC <- 4
    }
    #
    if(class(dev.list()) != "NULL"){dev.off()}
    #
    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    #
    for(j in 1:n_shocks){
        #
        if(save==TRUE){
            if(n_shocks==1){
                cairo_ps(filename="DSGEBIRFs.eps",height=height,width=width)
            }else{
                SaveIRF <- paste(varnames[nResp-n_shocks+j],"_Shock",".eps",sep="")
                #
                cairo_ps(filename=SaveIRF,height=height,width=width)
            }
        }
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))
        #
        IRFM <- IRFU <- IRFL <- Time <- NULL # CRAN check workaround
        VARIRFM <- VARIRFL <- VARIRFU <- NULL
        #
        if(n_shocks > 1 && observableIRFs!=TRUE){
            varnames2 <- varnames[-DropCount[-(j)]]
        }else{varnames2 <- varnames}
        #
        irf_plot_count <- 1
        IRFO <- 1
        if(observableIRFs==TRUE){IRFO <- nResp}
        #
        for(i in 1:MR){
            for(k in 1:MC){
                #
                if(irf_plot_count <= (nResp - n_shocks + IRFO)){
                    IRFDF <- data.frame(IRFDVPlot[,1:3,irf_plot_count,j],IRFDPlot[,,irf_plot_count,j])
                    colnames(IRFDF) <- c("VARIRFL","VARIRFM","VARIRFU","IRFL","IRFM","IRFU","Time")
                    #
                    if(comparison==TRUE){
                        print(ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(varnames2[irf_plot_count])) + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="darkslateblue",size=2) + geom_line(aes(y=VARIRFM),color="darkgreen",size=2) + geom_line(aes(y=VARIRFL),color="darkgreen",size=1,linetype=4) + geom_line(aes(y=VARIRFU),color="darkgreen",size=1,linetype=4) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                    }else{
                        print(ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(varnames2[irf_plot_count])) + geom_hline(yintercept=0) + geom_line(aes(y=VARIRFM),color="darkgreen",size=2) + geom_line(aes(y=VARIRFL),color="darkgreen",size=1,linetype=4) + geom_line(aes(y=VARIRFU),color="darkgreen",size=1,linetype=4) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                    }
                    #
                    irf_plot_count <- irf_plot_count + 1
                    #
                    Sys.sleep(0.3)
                }else{irf_plot_count <- irf_plot_count + 1}
            }
        }
        if(save==TRUE){dev.off()}
    }
    #
    return=list(IRFs=IRFsRet)
}
