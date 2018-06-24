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

IRF.Rcpp_bvarm <- function(obj,periods=10,cumulative=FALSE,cumul_inds=NULL,var_names=NULL,percentiles=c(.05,.50,.95),
                           which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,save=FALSE,save_format=c("pdf","eps"),
                           save_title=NULL,height=13,width=13,...)
{
    .irf_var(obj,periods,cumulative,cumul_inds,var_names,percentiles,which_shock,which_response,shocks_row_order,
             save,save_format,save_title,height,width)
}

IRF.Rcpp_bvars <- function(obj,periods=10,cumulative=FALSE,cumul_inds=NULL,var_names=NULL,percentiles=c(.05,.50,.95),
                           which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,save=FALSE,save_format=c("pdf","eps"),
                           save_title=NULL,height=13,width=13,...)
{
    .irf_var(obj,periods,cumulative,cumul_inds,var_names,percentiles,which_shock,which_response,shocks_row_order,
             save,save_format,save_title,height,width)
}

IRF.Rcpp_bvarcnw <- function(obj,periods=10,cumulative=FALSE,cumul_inds=NULL,var_names=NULL,percentiles=c(.05,.50,.95),
                             which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,save=FALSE,save_format=c("pdf","eps"),
                             save_title=NULL,height=13,width=13,...)
{
    .irf_var(obj,periods,cumulative,cumul_inds,var_names,percentiles,which_shock,which_response,shocks_row_order,
             save,save_format,save_title,height,width)
}

IRF.Rcpp_bvarinw <- function(obj,periods=10,cumulative=FALSE,cumul_inds=NULL,var_names=NULL,percentiles=c(.05,.50,.95),
                             which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,save=FALSE,save_format=c("pdf","eps"),
                             save_title=NULL,height=13,width=13,...)
{
    .irf_var(obj,periods,cumulative,cumul_inds,var_names,percentiles,which_shock,which_response,shocks_row_order,
             save,save_format,save_title,height,width)
}

IRF.Rcpp_cvar <- function(obj,periods=10,cumulative=FALSE,cumul_inds=NULL,var_names=NULL,percentiles=c(.05,.50,.95),
                          which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,save=FALSE,save_format=c("pdf","eps"),
                          save_title=NULL,height=13,width=13,...)
{
    .irf_var(obj,periods,cumulative,cumul_inds,var_names,percentiles,which_shock,which_response,shocks_row_order,
             save,save_format,save_title,height,width)
}

IRF.Rcpp_bvartvp <- function(obj,periods=10,which_irfs=NULL,var_names=NULL,percentiles=c(.05,.50,.95),
                             which_shock=NULL,which_response=NULL,save=FALSE,height=13,width=13,...)
{
    .irf_bvartvp(obj,periods,which_irfs,var_names,percentiles,which_shock,which_response,save,height,width)
}

IRF.Rcpp_gensys <- function(obj,periods=10,var_names=NULL,shocks_cov=NULL,save=FALSE,height=13,width=13,...)
{
    .irf_dsge(obj,periods,var_names,shocks_cov,save,height,width)
}

IRF.Rcpp_uhlig <- function(obj,periods=10,var_names=NULL,shocks_cov=NULL,save=FALSE,height=13,width=13,...)
{
    .irf_dsge(obj,periods,var_names,shocks_cov,save,height,width)
}

IRF.Rcpp_dsge_gensys <- function(obj,periods=10,obs_irfs=FALSE,var_names=NULL,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13,...)
{
    .irf_edsge(obj,periods,obs_irfs,var_names,percentiles,save,height,width)
}

IRF.Rcpp_dsge_uhlig <- function(obj,periods=10,obs_irfs=FALSE,var_names=NULL,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13,...)
{
    .irf_edsge(obj,periods,obs_irfs,var_names,percentiles,save,height,width)
}

IRF.Rcpp_dsgevar_gensys <- function(obj,periods,var_names=NULL,percentiles=c(.05,.50,.95),plot_comparison=TRUE,save=FALSE,height=13,width=13,...)
{
    .irf_dsgevar(obj,periods,var_names,percentiles,plot_comparison,save,height,width)
}

IRF.Rcpp_dsgevar_uhlig <- function(obj,periods,var_names=NULL,percentiles=c(.05,.50,.95),plot_comparison=TRUE,save=FALSE,height=13,width=13,...)
{
    .irf_dsgevar(obj,periods,var_names,percentiles,plot_comparison,save,height,width)
}

#

.irf_var <- function(obj, periods=10, cumulative=FALSE, cumul_inds=NULL, var_names=NULL, percentiles=c(.05,.50,.95), 
                     which_shock=NULL, which_response=NULL, shocks_row_order=TRUE, save=FALSE, save_format=c("pdf","eps"), 
                     save_title=NULL, height=13, width=13)
{
    #
    
    if (periods <= 0) {
        stop("error: need periods to be > 0")
    }

    #

    M <- obj$M
    n_draws <- dim(obj$beta_draws)[3]

    irf_temp <- obj$IRF(periods)$irf_vals

    # put the IRFs in a tesseract-type format

    irf_tess <- array(NA,dim=c(M,M,periods,n_draws))

    for (i in 1:n_draws) {
        irf_tess[,,,i] <- irf_temp[,,((i-1)*periods+1):(i*periods)]
    }
    
    if (cumulative)
    {
        if (is.null(cumul_inds)) {
            cumul_inds <- 1:M
        }

        for (i in 1:n_draws) {
            irfs_draw_i <- irf_tess[,,,i]
            
            for (ll in cumul_inds) {
                for (jj in 2:periods) {
                    irfs_draw_i[ll,,jj] <- irfs_draw_i[ll,,jj] + irfs_draw_i[ll,,jj-1]
                }
            }
            
            irf_tess[,,,i] <- irfs_draw_i
        }
    }

    irf_tess_1 <- apply(irf_tess,c(3,1,2),sort) # fix annoying bug
    irf_tess <- aperm(irf_tess_1,c(2,3,1,4))
    
    rm("irf_temp","irf_tess_1")

    irf_upper <- min(round(percentiles[3]*n_draws),n_draws)
    irf_mid <- round(percentiles[2]*n_draws)
    irf_lower <- max(round(percentiles[1]*n_draws),1)

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    
    #

    if (is.null(which_shock)) {
        which_shock <- 1:M
    }

    if (is.null(which_response)) {
        which_response <- 1:M
    }

    n_response <- length(which_response)
    n_shocks   <- length(which_shock)

    if (class(var_names) != "character") {
        var_names <- character(length=M)
        for (i in 1:M) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }

    #

    plot_vals <- array(NA,dim=c(periods,4,M,M))
    IRFPData <- 0

    for (i in 1:M) {
        for (k in 1:M) {
            IRFPData <- data.frame(irf_tess[,k,irf_lower,i],irf_tess[,k,irf_mid,i],irf_tess[,k,irf_upper,i],1:(periods))
            IRFPData <- as.matrix(IRFPData)
            plot_vals[,,k,i] <- IRFPData
        }
    }

    #
    # plot IRFs

    save_format <- match.arg(save_format)

    IRFL <- IRFM <- IRFU <- Time <- NULL # CRAN check workaround
    
    if (n_response == M && n_shocks == M) {

        if (save==TRUE)
        {
            if(class(dev.list()) != "NULL"){dev.off()}

            save_name <- ""
            if (!is.null(save_title)) {
                save_name <- paste(save_title,".",save_format,sep="")
            } else {
                save_name <- paste("IRFs.",save_format,sep="")
            }

            if (save_format=="eps") {
                cairo_ps(filename=save_name,height=height,width=width)
            } else {
                cairo_pdf(filename=save_name,height=height,width=width)
            }
        }
        
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(M,M)))
        
        for (i in 1:M) {
            for (k in 1:M) {
                NameResponse <- var_names[k]
                NameImpulse  <- var_names[i]
                
                IRFDF <- plot_vals[,,k,i]
                IRFDF <- data.frame(IRFDF)
                colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")
                
                #

                gg1 <- ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste("Shock from ",NameImpulse," to", NameResponse)) 
                gg2 <- gg1 + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="red",size=2) 
                gg3 <- gg2 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
                
                if (shocks_row_order) {
                    print(gg3,vp = vplayout(i,k))
                } else {
                    print(gg3,vp = vplayout(k,i))
                }
                
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
                if(class(dev.list()) != "NULL"){dev.off()}

                if (n_shocks==1) {
                    cairo_pdf(filename="IRFs.pdf",height=height,width=width)
                } else {
                    SaveIRF <- paste(var_names[i],"_Shock",".pdf",sep="")
                    cairo_pdf(filename=SaveIRF,height=height,width=width)
                }
            }

            if (save==TRUE)
            {
                if(class(dev.list()) != "NULL"){dev.off()}

                save_name <- ""
                if (n_shocks==1) {
                    if (!is.null(save_title)) {
                        save_name <- paste(save_title,".",save_format,sep="")
                    } else {
                        save_name <- paste("IRFs.",save_format,sep="")
                    }
                } else {
                    save_name <- paste(var_names[i],"_Shock",".",save_format,sep="")
                }

                if (save_format=="eps") {
                    cairo_ps(filename=save_name,height=height,width=width)
                } else {
                    cairo_pdf(filename=save_name,height=height,width=width)
                }
            }

            grid.newpage()
            pushViewport(viewport(layout=grid.layout(MR,MC)))
            
            for (k in which_response) {
                NameResponse <- var_names[k]
                NameImpulse  <- var_names[i]
                
                IRFDF <- plot_vals[,,k,i]
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

    return=list(plot_vals=plot_vals)
}

IRFcomp <- function(obj_1, obj_2, periods=10, cumulative=FALSE, cumul_inds=NULL, var_names=NULL, percentiles=c(.05,.50,.95), 
                    which_shock=NULL, which_response=NULL, shocks_row_order=TRUE, save=FALSE, save_format=c("pdf","eps"), 
                    save_title=NULL, height=13, width=13)
{
    #
    
    if (periods <= 0) {
        stop("error: need periods to be > 0")
    }

    #

    if (obj_1$M != obj_2$M) {
        stop("error: models need to be of the same data dimensions")
    }

    #

    M <- obj_1$M
    n_draws <- dim(obj_1$beta_draws)[3]

    irf_temp <- obj_1$IRF(periods)$irf_vals

    # put the IRFs in a tesseract-type format

    irf_tess <- array(NA,dim=c(M,M,periods,n_draws))

    for (i in 1:n_draws) {
        irf_tess[,,,i] <- irf_temp[,,((i-1)*periods+1):(i*periods)]
    }
    
    if (cumulative)
    {
        if (is.null(cumul_inds)) {
            cumul_inds <- 1:M
        }

        for (i in 1:n_draws) {
            irfs_draw_i <- irf_tess[,,,i]
            
            for (ll in cumul_inds) {
                for (jj in 2:periods) {
                    irfs_draw_i[ll,,jj] <- irfs_draw_i[ll,,jj] + irfs_draw_i[ll,,jj-1]
                }
            }
            
            irf_tess[,,,i] <- irfs_draw_i
        }
    }

    irf_tess_1 <- apply(irf_tess,c(3,1,2),sort) # fix annoying bug
    irf_tess <- aperm(irf_tess_1,c(2,3,1,4))
    
    rm("irf_temp","irf_tess_1")

    irf_upper <- min(round(percentiles[3]*n_draws),n_draws)
    irf_mid <- round(percentiles[2]*n_draws)
    irf_lower <- max(round(percentiles[1]*n_draws),1)

    #

    plot_vals_1 <- array(NA,dim=c(periods,4,M,M))
    IRFPData <- 0

    for (i in 1:M) {
        for (k in 1:M) {
            IRFPData <- data.frame(irf_tess[,k,irf_lower,i],irf_tess[,k,irf_mid,i],irf_tess[,k,irf_upper,i],1:(periods))
            IRFPData <- as.matrix(IRFPData)
            plot_vals_1[,,k,i] <- IRFPData
        }
    }

    #
    # Model 2

    M <- obj_2$M
    n_draws <- dim(obj_2$beta_draws)[3]

    irf_temp <- obj_2$IRF(periods)$irf_vals

    # put the IRFs in a tesseract-type format

    irf_tess <- array(NA,dim=c(M,M,periods,n_draws))

    for (i in 1:n_draws) {
        irf_tess[,,,i] <- irf_temp[,,((i-1)*periods+1):(i*periods)]
    }
    
    if (cumulative)
    {
        if (is.null(cumul_inds)) {
            cumul_inds <- 1:M
        }

        for (i in 1:n_draws) {
            irfs_draw_i <- irf_tess[,,,i]
            
            for (ll in cumul_inds) {
                for (jj in 2:periods) {
                    irfs_draw_i[ll,,jj] <- irfs_draw_i[ll,,jj] + irfs_draw_i[ll,,jj-1]
                }
            }
            
            irf_tess[,,,i] <- irfs_draw_i
        }
    }

    irf_tess_1 <- apply(irf_tess,c(3,1,2),sort) # fix annoying bug
    irf_tess <- aperm(irf_tess_1,c(2,3,1,4))
    
    rm("irf_temp","irf_tess_1")

    irf_upper <- min(round(percentiles[3]*n_draws),n_draws)
    irf_mid <- round(percentiles[2]*n_draws)
    irf_lower <- max(round(percentiles[1]*n_draws),1)

    #

    plot_vals_2 <- array(NA,dim=c(periods,4,M,M))
    IRFPData <- 0

    for (i in 1:M) {
        for (k in 1:M) {
            IRFPData <- data.frame(irf_tess[,k,irf_lower,i],irf_tess[,k,irf_mid,i],irf_tess[,k,irf_upper,i],1:(periods))
            IRFPData <- as.matrix(IRFPData)
            plot_vals_2[,,k,i] <- IRFPData
        }
    }

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    
    #

    if (is.null(which_shock)) {
        which_shock <- 1:M
    }

    if (is.null(which_response)) {
        which_response <- 1:M
    }

    n_response <- length(which_response)
    n_shocks   <- length(which_shock)

    if (class(var_names) != "character") {
        var_names <- character(length=M)
        for (i in 1:M) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }

    #
    # plot IRFs

    save_format <- match.arg(save_format)

    IRFL <- IRFM <- IRFU <- Time <- NULL # CRAN check workaround
    
    if (n_response == M && n_shocks == M) {

        if (save==TRUE)
        {
            if(class(dev.list()) != "NULL"){dev.off()}

            save_name <- ""
            if (!is.null(save_title)) {
                save_name <- paste(save_title,".",save_format,sep="")
            } else {
                save_name <- paste("IRFs.",save_format,sep="")
            }

            if (save_format=="eps") {
                cairo_ps(filename=save_name,height=height,width=width)
            } else {
                cairo_pdf(filename=save_name,height=height,width=width)
            }
        }
        
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(M,M)))
        
        for (i in 1:M) {
            for (k in 1:M) {
                NameResponse <- var_names[k]
                NameImpulse  <- var_names[i]
                
                IRFDF_1 <- plot_vals_1[,,k,i]
                #IRFDF_1 <- data.frame(IRFDF_1)

                IRFDF_2 <- plot_vals_2[,,k,i]
                IRFDF <- data.frame(IRFDF_1,IRFDF_2)

                #colnames(IRFDF_1) <- c("IRFL","IRFM","IRFU","Time")
                colnames(IRFDF) <- c("IRFL1","IRFM1","IRFU1","Time1","IRFL2","IRFM2","IRFU2","Time2")

                #

                gg1 <- ggplot(data=(IRFDF),aes(x=Time1)) + xlab("") + ylab(paste("Shock from ",NameImpulse," to", NameResponse)) + geom_hline(yintercept=0)
                gg2 <- gg1 + geom_line(aes(y=IRFM1),color="blue",size=2) + geom_line(aes(y=IRFL1),color="blue",size=1,linetype=4) + geom_line(aes(y=IRFU1),color="blue",size=1,linetype=4)
                gg3 <- gg2 + geom_line(aes(y=IRFM2),color="darkred",size=2) + geom_line(aes(y=IRFL2),color="darkred",size=1,linetype=4) + geom_line(aes(y=IRFU2),color="darkred",size=1,linetype=4)
                gg4 <- gg3 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))

                if (shocks_row_order) {
                    print(gg4,vp = vplayout(i,k))
                } else {
                    print(gg4,vp = vplayout(k,i))
                }
                
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
                if(class(dev.list()) != "NULL"){dev.off()}

                if (n_shocks==1) {
                    cairo_pdf(filename="IRFs.pdf",height=height,width=width)
                } else {
                    SaveIRF <- paste(var_names[i],"_Shock",".pdf",sep="")
                    cairo_pdf(filename=SaveIRF,height=height,width=width)
                }
            }

            if (save==TRUE)
            {
                if(class(dev.list()) != "NULL"){dev.off()}

                save_name <- ""
                if (n_shocks==1) {
                    if (!is.null(save_title)) {
                        save_name <- paste(save_title,".",save_format,sep="")
                    } else {
                        save_name <- paste("IRFs.",save_format,sep="")
                    }
                } else {
                    save_name <- paste(var_names[i],"_Shock",".",save_format,sep="")
                }

                if (save_format=="eps") {
                    cairo_ps(filename=save_name,height=height,width=width)
                } else {
                    cairo_pdf(filename=save_name,height=height,width=width)
                }
            }

            grid.newpage()
            pushViewport(viewport(layout=grid.layout(MR,MC)))
            
            for (k in which_response) {
                NameResponse <- var_names[k]
                NameImpulse  <- var_names[i]
                
                IRFDF_1 <- plot_vals_1[,,k,i]
                #IRFDF_1 <- data.frame(IRFDF_1)

                IRFDF_2 <- plot_vals_2[,,k,i]
                IRFDF <- data.frame(IRFDF_1,IRFDF_2)

                #colnames(IRFDF_1) <- c("IRFL","IRFM","IRFU","Time")
                colnames(IRFDF) <- c("IRFL1","IRFM1","IRFU1","Time1","IRFL2","IRFM2","IRFU2","Time2")

                #

                gg1 <- ggplot(data=(IRFDF),aes(x=Time1)) + xlab("") + ylab(paste("Shock from ", NameImpulse," to", NameResponse)) + geom_hline(yintercept=0)
                gg2 <- gg1 + geom_line(aes(y=IRFM1),color="blue",size=2) + geom_line(aes(y=IRFL1),color="blue",size=1,linetype=4) + geom_line(aes(y=IRFU1),color="blue",size=1,linetype=4)
                gg3 <- gg2 + geom_line(aes(y=IRFM2),color="darkred",size=2) + geom_line(aes(y=IRFL2),color="darkred",size=1,linetype=4) + geom_line(aes(y=IRFU2),color="darkred",size=1,linetype=4)
                gg4 <- gg3 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))

                print(gg4,vp = vplayout(plot_ind_r,plot_ind_c))
                
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

    return=list(plot_vals_1=plot_vals_1,plot_vals_2=plot_vals_2)
}

.irf_bvartvp <- function(obj, periods=10, which_irfs=NULL, var_names=NULL, percentiles=c(.05,.50,.95), 
                         which_shock=NULL, which_response=NULL, save=FALSE, height=13, width=13)
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
    
    #

    if (is.null(which_shock)) {
        which_shock <- 1:M
    }

    if (is.null(which_response)) {
        which_response <- 1:M
    }

    n_response <- length(which_response)
    n_shocks   <- length(which_shock)

    if (class(var_names) != "character") {
        var_names <- character(length=M)
        for (i in 1:M) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }

    #
    # loop

    if (is.null(which_irfs)) {
        which_irfs <- dim(obj$alpha_draws)[2]
    }

    n_plots <- length(which_irfs)

    for (pp in 1:n_plots) {

        irf_temp <- obj$IRF(periods,which_irfs[pp]-1)$irf_vals # minus 1 to account for zero-index

        # put the IRFs in a tesseract-type format
        
        irf_tess <- array(NA,dim=c(M,M,periods,n_draws))

        for (i in 1:n_draws) {
            irf_tess[,,,i] <- irf_temp[,,((i-1)*periods+1):(i*periods)]
        }

        rm("irf_temp")

        irf_tess <- apply(irf_tess,c(3,1,2),sort)
        irf_tess <- aperm(irf_tess,c(2,3,1,4))
        
        plot_vals <- array(NA,dim=c(periods,4,M,M))

        for (i in 1:M) {
            for (k in 1:M) {
                IRFPData <- data.frame(irf_tess[,k,irf_lower,i],irf_tess[,k,irf_mid,i],irf_tess[,k,irf_upper,i],1:(periods))
                IRFPData <- as.matrix(IRFPData)
                plot_vals[,,k,i] <- IRFPData
            }
        }

        #
        # plot IRFs

        IRFL <- IRFM <- IRFU <- Time <- NULL # CRAN check workaround
        
        if (n_response == M && n_shocks == M) {

            if (save==TRUE) {
                if(class(dev.list()) != "NULL"){dev.off()}
                cairo_ps(filename="IRFs.eps",height=height,width=width)
            }
            
            grid.newpage()
            pushViewport(viewport(layout=grid.layout(M,M)))
            
            for (i in 1:M) {
                for (k in 1:M) {
                    NameResponse <- var_names[k]
                    NameImpulse  <- var_names[i]
                    
                    IRFDF <- plot_vals[,,k,i]
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
                    if(class(dev.list()) != "NULL"){dev.off()}

                    if (n_shocks==1) {
                        cairo_ps(filename="IRFs.eps",height=height,width=width)
                    } else {
                        SaveIRF <- paste(var_names[i],"_Shock",".eps",sep="")
                        cairo_ps(filename=SaveIRF,height=height,width=width)
                    }
                }

                grid.newpage()
                pushViewport(viewport(layout=grid.layout(MR,MC)))
                
                for (k in which_response) {
                    NameResponse <- var_names[k]
                    NameImpulse  <- var_names[i]
                    
                    IRFDF <- plot_vals[,,k,i]
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

    return=list(plot_vals=plot_vals)
}

#
# DSGE

.irf_dsge <- function(obj,periods=10,var_names=NULL,shocks_cov=NULL,save=FALSE,height=13,width=13)
{
    
    if (periods <= 0) {
        stop("error: need periods > 0")
    }

    #

    irfs <- obj$IRF(periods)$irf_vals

    irfs <- round(irfs,10)

    M <- dim(irfs)[2]
    n_shocks <- dim(irfs)[3]

    if (class(var_names) != "character") {
        var_names <- character(length=M)
        for (i in 1:M) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }
    
    #
    
    plot_vals <- irfs
    drop_ind <- (M - n_shocks + 1):M
    
    if (n_shocks > 1) { # drop other shocks
        plot_vals <- array(0, dim=c(periods, M - n_shocks + 1, n_shocks))
        
        for (j in 1:n_shocks) {
            plot_vals[,,j] <- irfs[,-drop_ind[-(j)],j]
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
    
    #

    var_names2 <- var_names

    for (j in 1:n_shocks) {
        
        if (save==TRUE) {
            if(class(dev.list()) != "NULL"){dev.off()}

            if (n_shocks==1) {
                cairo_ps(filename="DSGEIRFs.eps",height=height,width=width)
            } else {
                SaveIRF <- paste(var_names[M-n_shocks+j],"_Shock",".eps",sep="")
                cairo_ps(filename=SaveIRF,height=height,width=width)
            }
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))

        if (n_shocks > 1) {
            var_names2 <- var_names[-drop_ind[-(j)]]
        } else {
            var_names2 <- var_names
        }

        #

        plot_count <- 1

        for (i in 1:MR) {
            for (k in 1:MC) {
                
                if (plot_count <= (M - n_shocks + 1)) {
                    IRFDF <- data.frame(plot_vals[,plot_count,j],1:periods)
                    colnames(IRFDF) <- c("IRFM","Time")
                    
                    print(ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(var_names2[plot_count])) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="darkslateblue",size=2) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                    
                    #

                    plot_count <- plot_count + 1
                    
                    Sys.sleep(0.3)
                } else {
                    plot_count <- plot_count + 1
                }
            }
        }

        if(save==TRUE){dev.off()}
    }
    #
    return=list(plot_vals=plot_vals)
}

.irf_edsge <- function(obj,periods=10,obs_irfs=FALSE,var_names=NULL,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13)
{    
    if (periods <= 0) {
        stop("error: need periods > 0")
    }

    #

    n_draws <- dim(obj$dsge_draws)[1]

    if (n_draws <= 0) {
        stop("error: no MCMC draws detected")
    }

    irfs <- obj$IRF(periods,obs_irfs)$irf_vals # deep copy is needed
    irfs <- round(irfs,10)

    M <- dim(irfs)[2]
    n_shocks <- dim(irfs)[3] / n_draws

    n_response <- M - n_shocks + 1
    
    if (obs_irfs) {
        n_response <- M
    }

    #

    irfs_i <- array(0, dim=c(periods, M, n_shocks))
    irfs_tess <- array(0, dim=c(periods, n_response, n_shocks, n_draws))

    drop_ind <- n_response:M
    
    for (i in 1:n_draws) {
        irfs_i <- irfs[,,((i-1)*n_shocks + 1):(i*n_shocks)]

        if (n_shocks > 1 && obs_irfs==FALSE) {
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

    plot_vals <- array(NA,dim=c(periods,4,n_response,n_shocks))

    for (i in 1:n_shocks) {
        for (k in 1:n_response) {
            plot_vals_temp <- data.frame(c(irfs_tess[irf_lower,,k,i]),c(irfs_tess[irf_mid,,k,i]),c(irfs_tess[irf_upper,,k,i]),1:(periods))
            plot_vals[,,k,i] <- data.matrix(plot_vals_temp)
        }
    }

    #

    if (class(var_names) != "character") {
        var_names <- character(length=M)
        for (i in 1:M) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }

    #

    IRFM <- IRFU <- IRFL <- Time <- NULL # CRAN check workaround

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
    
    #

    var_names2 <- var_names

    for (j in 1:n_shocks) {
        #
        if (save==TRUE) {
            if(class(dev.list()) != "NULL"){dev.off()}

            if (n_shocks==1) {
                cairo_ps(filename="DSGEBIRFs.eps",height=height,width=width)
            } else {
                SaveIRF <- paste(var_names[M-n_shocks+j],"_Shock",".eps",sep="")
                #
                cairo_ps(filename=SaveIRF,height=height,width=width)
            }
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))

        #

        if (n_shocks > 1 && obs_irfs==FALSE) {
            var_names2 <- var_names[-drop_ind[-(j)]]
        }

        #

        plot_count <- 1

        for (i in 1:MR) {
            for (k in 1:MC) {
                #
                if (plot_count <= n_response) {
                    IRFDF <- data.frame(plot_vals[,,plot_count,j])
                    colnames(IRFDF) <- c("IRFL","IRFM","IRFU","Time")

                    #

                    gg1 <- ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(var_names2[plot_count])) 
                    gg2 <- gg1 + geom_ribbon(aes(ymin=IRFL,ymax=IRFU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=IRFM),color="darkslateblue",size=2) 
                    gg3 <- gg2 + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89'))
                    print(gg3,vp = vplayout(i,k))
                    
                    #

                    plot_count <- plot_count + 1
                    
                    Sys.sleep(0.3)
                } else {
                    plot_count <- plot_count + 1
                }
            }
        }

        if(save==TRUE){dev.off()}
    }
    #
    return=list(plot_vals=plot_vals)
}

.irf_dsgevar <- function(obj,periods,var_names=NULL,percentiles=c(.05,.50,.95),comparison_plot=TRUE,save=FALSE,height=13,width=13)
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

    #

    irf_upper <- round(percentiles[3]*n_draws)
    irf_mid <- round(percentiles[2]*n_draws)
    irf_lower <- round(percentiles[1]*n_draws)

    #

    plot_vals <- array(NA,dim=c(periods,4,n_response,n_shocks))

    for (i in 1:M) {
        for (k in 1:M) {
            plot_vals_temp <- data.frame(c(irfs_tess[irf_lower,,k,i]),c(irfs_tess[irf_mid,,k,i]),c(irfs_tess[irf_upper,,k,i]),1:(periods))
            plot_vals[,,k,i] <- data.matrix(plot_vals_temp)
        }
    }
    
    #

    if (class(var_names) != "character") {
        var_names <- character(length=M)
        for (i in 1:M) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }

    #

    IRFM <- IRFU <- IRFL <- Time <- NULL # CRAN check workaround

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
    
    #

    for (j in 1:n_shocks) {
        #
        if (save==TRUE) {
            if(class(dev.list()) != "NULL"){dev.off()}

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

        plot_count <- 1

        for (i in 1:MR) {
            for (k in 1:MC) {
                #
                if (plot_count <= n_response) {
                    
                    IRFDF <- data.frame(plot_vals[,,plot_count,j])
                    colnames(IRFDF) <- c("VARIRFL","VARIRFM","VARIRFU","Time")

                    #
                    
                    print(ggplot(data=(IRFDF),aes(x=Time)) + xlab("") + ylab(paste(var_names[plot_count])) + geom_hline(yintercept=0) + geom_line(aes(y=VARIRFM),color="darkgreen",size=2) + geom_line(aes(y=VARIRFL),color="darkgreen",size=1,linetype=4) + geom_line(aes(y=VARIRFU),color="darkgreen",size=1,linetype=4) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))

                    #
                    
                    plot_count <- plot_count + 1
                    Sys.sleep(0.3)
                    
                } else {
                    plot_count <- plot_count + 1
                }
            }
        }

        if(save==TRUE){dev.off()}

    }

    #

    return=list(plot_vals=plot_vals)
}
