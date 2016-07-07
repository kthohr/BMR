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

gacf<-function(y,lags=10,ci=.95,plot=TRUE,barcolor="purple",names=FALSE,save=FALSE,height=12,width=12){
    .ggplotacf(y,lags,ci,plot,barcolor,names,save,height,width)
}

gpacf<-function(y,lags=10,ci=.95,plot=TRUE,barcolor="darkred",names=FALSE,save=FALSE,height=12,width=12){
    .ggplotpacf(y,lags,ci,plot,barcolor,names,save,height,width)
}

.ggplotacf<-function(y,lags=10,ci=.95,plot=TRUE,barcolor="purple",names=FALSE,save=FALSE,height=12,width=12){
    #
    M <- as.numeric(ncol(y))
    #
    CIMat <- matrix(NA,lags,M)
    ACFMat <- matrix(NA,lags,M)
    #
    if(class(colnames(y)) != "character"){
        varnames <- ncol(y)
        for(i in 1:ncol(y)){  
            varnames[i] <- paste("Var",i,sep="")
        }
    }else{
        varnames <- colnames(y)
    }
    #
    for(j in 1:M){
        myacf <- acf(y[,j],lags,plot=FALSE)
        myacf2 <- as.numeric(myacf$acf)
        myacf2 <- myacf2[2:(lags+1)]
        ACFMat[,j] <- myacf2
        myci <- qnorm((1 + ci)/2)/sqrt(sum(!is.na(y[,j])))
        mycount <- 1:lags
        CIMat[1,j] <- myci
        for(i in 2:(length(mycount))){
            CIMat[i,j] <- qnorm((1 + ci)/2)*sqrt((1+2*sum(t(myacf2[1:i-1])*myacf2[1:i-1]))/sum(!is.na(y[,j])))
        }
    }
    #
    vplayout<-function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    #
    if(plot==TRUE){
        if(save==TRUE){
            if(class(dev.list()) != "NULL"){dev.off()}
            cairo_ps(filename="ACF.eps",height=height,width=width)
        }else{
            grid.newpage()
        }
        pushViewport(viewport(layout=grid.layout(1,M)))
        for(j in 1:M){
            kacf <- ACFMat[,j]
            kci <- CIMat[,j]
            VarName <- varnames[j]
            #
            ACFdata <- data.frame(mycount,kacf,kci)
            #
            aa <- ggplot(ACFdata, aes(x=mycount,y=kacf))
            aa <- aa + geom_bar(fill = I(paste(barcolor)),stat="identity")
            #
            if(sum(kacf < 0) == 0){
                if(names==TRUE){
                    aaa <- aa + geom_line(aes(y=kci),lty=2,color="red") + geom_ribbon(aes(ymin=0,ymax=kci),fill="blue",alpha=0.15) + xlab("Lags") + ylab("") + labs(title=VarName)
                }else if(names==FALSE){
                    aaa <- aa + geom_line(aes(y=kci),lty=2,color="red") + geom_ribbon(aes(ymin=0,ymax=kci),fill="blue",alpha=0.15) + xlab("Lags") + ylab("")
                }
            }else{
                if(names==TRUE){
                    aaa <- aa + geom_line(aes(y=kci),lty=2,color="red") + geom_line(aes(y=-kci),lty=2,color="red") + geom_ribbon(aes(ymin=-kci,ymax=kci),fill="blue",alpha=0.15) + xlab("Lags") + ylab("") + labs(title=VarName)
                }else if(names==FALSE){
                    aaa <- aa + geom_line(aes(y=kci),lty=2,color="red") + geom_line(aes(y=-kci),lty=2,color="red") + geom_ribbon(aes(ymin=-kci,ymax=kci),fill="blue",alpha=0.15) + xlab("Lags") + ylab("")
                }
            }
            #
            suppressWarnings(print(aaa,vp = vplayout(1,j)))
        }
        if(save==TRUE){dev.off()}
    }
    #
}

.ggplotpacf<-function(y,lags=10,ci=.95,plot=TRUE,barcolor="darkred",names=FALSE,save=FALSE,height=12,width=12){
    #
    M <- as.numeric(ncol(y))
    #
    CIMat <- matrix(NA,lags,M)
    PACFMat <- matrix(NA,lags,M)
    #
    if(class(colnames(y)) != "character"){
        varnames <- ncol(y)
        for(i in 1:ncol(y)){  
            varnames[i] <- paste("Var",i,sep="")
        }
    }else{
        varnames <- colnames(y)
    }
    #
    for(j in 1:M){
        mypacf <- pacf(y[,j],lags,plot=FALSE)
        mypacf2 <- as.numeric(mypacf$acf)
        mypacf2 <- mypacf2[1:lags]
        PACFMat[,j] <- mypacf2
        myci <- qnorm((1 + ci)/2)/sqrt(sum(!is.na(y[,j])))
        mycount <- 1:lags
        CIMat[1,j] <- myci
        for(i in 2:(length(mycount))){
            CIMat[i,j] <- qnorm((1 + ci)/2)*sqrt((1+2*sum(t(mypacf2[1:i-1])*mypacf2[1:i-1]))/sum(!is.na(y[,j])))
        }
    }
    #
    vplayout<-function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    #
    if(plot==TRUE){
        if(save==TRUE){
            if(class(dev.list()) != "NULL"){dev.off()}
            cairo_ps(filename="PACF.eps",height=height,width=width)
        }else{
            grid.newpage()
        }
        pushViewport(viewport(layout=grid.layout(1,M)))
        for(j in 1:M){
            kpacf <- PACFMat[,j]
            kci <- CIMat[,j]
            VarName <- varnames[j]
            #
            PACFdata <- data.frame(mycount,kpacf,kci)
            #
            aa <- ggplot(PACFdata, aes(x=mycount,y=kpacf))
            aa <- aa + geom_bar(fill = I(paste(barcolor)),stat="identity")
            #
            if(sum(kpacf < 0) == 0){
                if(names==TRUE){
                    aaa <- aa + geom_line(aes(y=kci),lty=2,color="red") + geom_ribbon(aes(ymin=0,ymax=kci),fill="blue",alpha=0.15) + xlab("Lags") + ylab("") + labs(title=VarName)
                }else if(names==FALSE){
                    aaa <- aa + geom_line(aes(y=kci),lty=2,color="red") + geom_ribbon(aes(ymin=0,ymax=kci),fill="blue",alpha=0.15) + xlab("Lags") + ylab("")
                }
            }else{
                if(names==TRUE){
                    aaa <- aa + geom_line(aes(y=kci),lty=2,color="red") + geom_line(aes(y=-kci),lty=2,color="red") + geom_ribbon(aes(ymin=-kci,ymax=kci),fill="blue",alpha=0.15) + xlab("Lags") + ylab("") + labs(title=VarName)
                }else if(names==FALSE){
                    aaa <- aa + geom_line(aes(y=kci),lty=2,color="red") + geom_line(aes(y=-kci),lty=2,color="red") + geom_ribbon(aes(ymin=-kci,ymax=kci),fill="blue",alpha=0.15) + xlab("Lags") + ylab("")
                }
            }
            #
            suppressWarnings(print(aaa,vp = vplayout(1,j)))
        }
        if(save==TRUE){dev.off()}
    }
    #
}