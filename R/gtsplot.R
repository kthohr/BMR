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

gtsplot <- function(X,dates=NULL,rowdates=FALSE,dates.format="%Y-%m-%d",save=FALSE,height=13,width=11){
    .bmrgtsplot(X,dates,rowdates,dates.format,save,height,width)
}

.bmrgtsplot <- function(X,dates=NULL,rowdates=FALSE,dates.format="%Y-%m-%d",save=FALSE,height=13,width=11){
    #
    nSeries <- ncol(X)
    if(class(nSeries) == "NULL"){
        nSeries <- 1
    }else{
        nSeries <- as.numeric(nSeries)
    }
    #
    if(rowdates==T){
        XDates <- rownames(X)
    }
    if(class(dates)!="NULL"){
        XDates <- dates
    }
    #
    XMat <- as.matrix(X)
    #
    if(rowdates==T || class(dates)!="NULL"){
        dates <- as.Date(XDates, dates.format)
    }else{
        dates <- 1:nrow(XMat)
    }
    #
    if(class(colnames(XMat)) != "character"){
        VarNames <- character(length=nSeries)
        for(i in 1:nSeries){  
            VarNames[i] <- paste("Variable",i,sep="")
        }
    }else{
        VarNames <-colnames(X)
    }
    #
    MR <- 0; MC <- 0
    plotpages <- 1
    #
    if(nSeries < 4){
        MR <- nSeries; MC <- 1
    }else if(nSeries == 4){
        MR <- 2; MC <-2
    }else if(nSeries > 4 && nSeries < 7){
        MR <- 3; MC <- 2
    }else if(nSeries > 6 && nSeries < 10){
        MR <- 3; MC <- 3
    }else if(nSeries > 9 && nSeries < 13){
        MR <- 4; MC <- 3
    }else if(nSeries > 12 && nSeries < 25){
        MR <- 4; MC <- 3
        plotpages <- 2
    }else if(nSeries > 24 && nSeries < 37){
        MR <- 4; MC <- 3
        plotpages <- 3
    }else if(nSeries > 36 && nSeries < 49){
        MR <- 4; MC <- 3
        plotpages <- 4
    }else if(nSeries > 48 && nSeries < 61){
        MR <- 4; MC <- 3
        plotpages <- 5
    }else if(nSeries > 60 && nSeries < 73){
        MR <- 4; MC <- 3
        plotpages <- 6
    }
    #
    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    #
    SeriesCount <- 1
    Orig <- NULL # CRAN check workaround
    #
    for(j in 1:plotpages){
        #
        if(save==TRUE){
            if(class(dev.list()) != "NULL"){dev.off()}
            if(plotpages==1){
                cairo_ps(filename="TimeSeriesPlot.eps",height=height,width=width)
            }else{
                SaveIRF <- paste("TimeSeries_",j,".eps",sep="")
                #
                cairo_ps(filename=SaveIRF,height=height,width=width)
            }
        }else{
            grid.newpage()
        }
        pushViewport(viewport(layout=grid.layout(MR,MC)))
        #
        for(i in 1:MR){
            for(k in 1:MC){
                #
                if(SeriesCount <= nSeries){
                    GGMat <- data.frame(dates,XMat[,SeriesCount])
                    colnames(GGMat) <- c("dates","Orig")
                    VarName <- VarNames[SeriesCount]
                    #
                    print(ggplot(GGMat,aes(dates,Orig)) + xlab(NULL) + ylab(paste(VarName)) + geom_line(colour="royalblue4") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                    #
                    SeriesCount <- SeriesCount + 1
                    #
                    Sys.sleep(0.6)
                }else{SeriesCount <- SeriesCount + 1}
            }
        }
        if(save==TRUE){dev.off()}
    }
    #
}