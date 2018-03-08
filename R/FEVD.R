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

FEVD.Rcpp_bvarm <- function(obj,periods=10,var_names=NULL,percentiles=c(.05,.50,.95),
                            which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,
                            save=FALSE,save_format=c("pdf","eps"),
                            save_title=NULL,height=13,width=13,...)
{
    .fevd_var(obj,periods,var_names,percentiles,which_shock,which_response,shocks_row_order,
              save,save_format,save_title,height,width)
}

FEVD.Rcpp_bvars <- function(obj,periods=10,var_names=NULL,percentiles=c(.05,.50,.95),
                            which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,
                            save=FALSE,save_format=c("pdf","eps"),
                            save_title=NULL,height=13,width=13,...)
{
    .fevd_var(obj,periods,var_names,percentiles,which_shock,which_response,shocks_row_order,
              save,save_format,save_title,height,width)
}

FEVD.Rcpp_bvarcnw <- function(obj,periods=10,var_names=NULL,percentiles=c(.05,.50,.95),
                              which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,
                              save=FALSE,save_format=c("pdf","eps"),
                              save_title=NULL,height=13,width=13,...)
{
    .fevd_var(obj,periods,var_names,percentiles,which_shock,which_response,shocks_row_order,
              save,save_format,save_title,height,width)
}

FEVD.Rcpp_bvarinw <- function(obj,periods=10,var_names=NULL,percentiles=c(.05,.50,.95),
                              which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,
                              save=FALSE,save_format=c("pdf","eps"),
                              save_title=NULL,height=13,width=13,...)
{
    .fevd_var(obj,periods,var_names,percentiles,which_shock,which_response,shocks_row_order,
              save,save_format,save_title,height,width)
}

FEVD.Rcpp_cvar <- function(obj,periods=10,var_names=NULL,percentiles=c(.05,.50,.95),
                           which_shock=NULL,which_response=NULL,shocks_row_order=TRUE,
                           save=FALSE,save_format=c("pdf","eps"),
                           save_title=NULL,height=13,width=13,...)
{
    .fevd_var(obj,periods,var_names,percentiles,which_shock,which_response,shocks_row_order,
              save,save_format,save_title,height,width)
}

#

.fevd_var <- function(obj, periods=10, var_names=NULL, percentiles=c(.05,.50,.95), 
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

    fevd_temp <- obj$FEVD(periods)$fevd_vals

    # put the FEVDs in a tesseract-type format

    fevd_tess <- array(NA,dim=c(M,M,periods,n_draws))

    for (i in 1:n_draws) {
        fevd_tess[,,,i] <- fevd_temp[,,((i-1)*periods+1):(i*periods)]
    }

    rm("fevd_temp")

    fevd_tess <- apply(fevd_tess,c(3,1,2),sort)
    fevd_tess <- aperm(fevd_tess,c(2,3,1,4))

    fevd_upper <- round(percentiles[3]*n_draws)
    fevd_mid <- round(percentiles[2]*n_draws)
    fevd_lower <- round(percentiles[1]*n_draws)

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
    FEVDPData <- 0

    for (i in 1:M) {
        for (k in 1:M) {
            FEVDPData <- data.frame(fevd_tess[,k,fevd_lower,i],fevd_tess[,k,fevd_mid,i],fevd_tess[,k,fevd_upper,i],1:(periods))
            FEVDPData <- as.matrix(FEVDPData)
            plot_vals[,,k,i] <- FEVDPData
        }
    }

    #
    # plot FEVDs

    save_format <- match.arg(save_format)

    FEVDL <- FEVDM <- FEVDU <- Time <- NULL # CRAN check workaround
    
    if (n_response == M && n_shocks == M) {

        if (save==TRUE)
        {
            if(class(dev.list()) != "NULL"){dev.off()}

            save_name <- ""
            if (!is.null(save_title)) {
                save_name <- paste(save_title,".",save_format,sep="")
            } else {
                save_name <- paste("FEVDs.",save_format,sep="")
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
                
                FEVDDF <- plot_vals[,,k,i]
                FEVDDF <- data.frame(FEVDDF)
                colnames(FEVDDF) <- c("FEVDL","FEVDM","FEVDU","Time")
                
                #

                gg1 <- ggplot(data=(FEVDDF),aes(x=Time)) + xlab("") + ylab(paste("Shock from ",NameImpulse," to", NameResponse)) 
                gg2 <- gg1 + geom_ribbon(aes(ymin=FEVDL,ymax=FEVDU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=FEVDM),color="red",size=2) 
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
            stop("You have too many FEVDs to plot!")
        }

        #
        
        for (i in which_shock) {
            plot_ind_r <- 1
            plot_ind_c <- 1
            
            if (save==TRUE) {
                if(class(dev.list()) != "NULL"){dev.off()}

                if (n_shocks==1) {
                    cairo_pdf(filename="FEVDs.pdf",height=height,width=width)
                } else {
                    SaveFEVD <- paste(var_names[i],"_Shock",".pdf",sep="")
                    cairo_pdf(filename=SaveFEVD,height=height,width=width)
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
                        save_name <- paste("FEVDs.",save_format,sep="")
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
                
                FEVDDF <- plot_vals[,,k,i]
                FEVDDF <- data.frame(FEVDDF)

                colnames(FEVDDF) <- c("FEVDL","FEVDM","FEVDU","Time")

                #

                gg1 <- ggplot(data=(FEVDDF),aes(x=Time)) + xlab("") + ylab(paste("Shock from ", NameImpulse," to", NameResponse)) 
                gg2 <- gg1 + geom_ribbon(aes(ymin=FEVDL,ymax=FEVDU),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_hline(yintercept=0) + geom_line(aes(y=FEVDM),color="red",size=2) 
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
