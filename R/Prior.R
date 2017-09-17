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

prior <- function(priorform,priorpars,parname=NULL,moments=TRUE,NR=NULL,NC=NULL)
{
    .prior_function(priorform,priorpars,parname,moments,NR,NC)
}

.prior_function <- function(priorform,priorpars,parname=NULL,moments=TRUE,NR=NULL,NC=NULL)
{
    #
    x <- 0; PriorDF <- 0;
    xmean <- 0; xmode <- 0; xvar <- 0
    #
    if (priorform=="beta") {
        #
        x <- seq(0, 1, length=400)
        betaA <- priorpars[1]; betaB <- priorpars[2]
        hx <- dbeta(x,betaA,betaB)
        #
        xmean <- (betaA)/(betaA + betaB)
        xmode <- NaN
        if(betaA > 1 && betaB > 1){
            xmode <- (betaA - 1)/(betaA + betaB - 2)
        }
        xvar <- (betaA*betaB)/( ((betaA + betaB)^2)*(betaA + betaB + 1) )
        #
        if (moments==TRUE) {
            cat('Mean: ', xmean, '. \n', sep="")
            cat('Mode: ', xmode, '. \n', sep="")
            cat('Variance: ', xvar, '. \n', sep="")
        }
        #
        PriorDF <- data.frame(hx,x)
        colnames(PriorDF) <- c("Density","Domain")
        #
    } else if (priorform=="uniform") {
        #
        par1 <- priorpars[1]; par2 <- priorpars[2]
        x <- seq(par1 - (par2 - par1)/10, par2 + (par2 - par1)/10, length=400)
        hx <- dunif(x,par1,par2)
        #
        xmean <- 0.5*(par1 + par2)
        xmode <- paste("any value in the continuum between ", par1, " and ", par2,sep="")
        xvar <- (1/12)*(par2 - par1)^2
        #
        if(moments==TRUE){
            cat('Mean: ', xmean, '. \n', sep="")
            cat('Mode: ', xmode, '. \n', sep="")
            cat('Variance: ', xvar, '. \n', sep="")
        }
        #
        PriorDF <- data.frame(hx,x)
        colnames(PriorDF) <- c("Density","Domain")
        #
    } else if (priorform=="igamma") {
        #
        InvGamma <- function(x,par1,par2){
            dens <- (par2^(par1)/gamma(par1))*(x^(-par1-1))*exp(-par2/x)
            return(dens)
        }
        #
        x <- seq(0.01, 5, length=400)
        par1 <- priorpars[1]; par2 <- priorpars[2]
        hx <- InvGamma(x,par1,par2)
        #
        xmean <- (par2)/(par1 - 1)
        xmode <- (par2)/(par1 + 1)
        xvar <- NaN
        if(par1 > 2){
            xvar <- (par2^2)/( ((par1 - 1)^2)*(par1 - 2) )
        }
        #
        if(moments==TRUE){
            cat('Mean: ', xmean, '. \n', sep="")
            cat('Mode: ', xmode, '. \n', sep="")
            cat('Variance: ', xvar, '. \n', sep="")
        }
        #
        PriorDF <- data.frame(hx,x)
        colnames(PriorDF) <- c("Density","Domain")
        #
    } else if (priorform=="gamma") {
        #
        par1 <- priorpars[1]; par2 <- 1/priorpars[2]
        if(par2 < 1){
            x <- seq(0.01, 2*par1/(par2^2), length=400)
        }else{
            x <- seq(0.01, par1/par2 + 2*par2*par1/(par2^2), length=400)
        }
        hx <- dgamma(x,shape=par1,rate=par2)
        #
        xmean <- par1/par2
        xmode <- (par1 - 1)/par2
        xvar <- par1/(par2^2)
        #
        if(moments==TRUE){
            cat('Mean: ', xmean, '. \n', sep="")
            cat('Mode: ', xmode, '. \n', sep="")
            cat('Variance: ', xvar, '. \n', sep="")
        }
        #
        PriorDF <- data.frame(hx,x)
        colnames(PriorDF) <- c("Density","Domain")
        #
    } else if(priorform=="normal") {
        #
        par1 <- priorpars[1]; par2 <- sqrt(priorpars[2])
        x <- seq(par1-par2*3, par1+par2*3, length=400)
        hx <- dnorm(x,par1,par2)
        #
        xmean <- par1
        xmode <- par1
        xvar <- par2^2
        #
        if(moments==TRUE){
            cat('Mean: ', xmean, '. \n', sep="")
            cat('Mode: ', xmode, '. \n', sep="")
            cat('Variance: ', xvar, '. \n', sep="")
        }
        #
        PriorDF <- data.frame(hx,x)
        colnames(PriorDF) <- c("Density","Domain")
        #
    } else {
        stop("unknown prior form.\n",call.=FALSE)
    }
    
    #
    #

    Domain <- Density <- NULL # CRAN check workaround
    
    if(class(NR)=="NULL" | class(NC)=="NULL"){
        if(class(parname)=="character"){
            if(priorform=="uniform"){
                print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_line(aes(y=Density),color="red",size=1) + labs(title=parname))
            }else{
                print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_vline(xintercept=xmode,linetype = "longdash") + geom_line(aes(y=Density),color="red",size=1) + labs(title=parname))
            }
        }else{
            if(priorform=="uniform"){
                print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_line(aes(y=Density),color="red",size=1))
            }else{
                print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_vline(xintercept=xmode,linetype = "longdash") + geom_line(aes(y=Density),color="red",size=1))
            }
        }
    }else{
        #
        vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
        #
        if(class(parname)=="character"){
            if(priorform=="uniform"){
                print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_line(aes(y=Density),color="red",size=1) + labs(title=parname),vp = vplayout(NR,NC))
            }else{
                print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_vline(xintercept=xmode,linetype = "longdash") + geom_line(aes(y=Density),color="red",size=1) + labs(title=parname),vp = vplayout(NR,NC))
            }
        }else{
            if(priorform=="uniform"){
                print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_line(aes(y=Density),color="red",size=1),vp = vplayout(NR,NC))
            }else{
                print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_vline(xintercept=xmode,linetype = "longdash") + geom_line(aes(y=Density),color="red",size=1),vp = vplayout(NR,NC))
            }
        }
    }
}
