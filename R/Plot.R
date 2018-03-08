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

plot.Rcpp_bvarm <- function(x,type=1,var_names=NULL,save=FALSE,height=13,width=13,...)
{
    .plotbvar(x,type,var_names,save,height,width)
}

plot.Rcpp_bvars <- function(x,type=1,var_names=NULL,save=FALSE,height=13,width=13,...)
{
    .plotbvar(x,type,var_names,save,height,width)
}

plot.Rcpp_bvarcnw <- function(x,type=1,var_names=NULL,save=FALSE,height=13,width=13,...)
{
    .plotbvar(x,type,var_names,save,height,width)
}

plot.Rcpp_bvarinw <- function(x,type=1,var_names=NULL,save=FALSE,height=13,width=13,...)
{
    .plotbvar(x,type,var_names,save,height,width)
}

plot.Rcpp_bvartvp <- function(x,var_names=NULL,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13,...)
{
    .plotbvartvp(x,var_names,percentiles,save,height,width)
}

plot.Rcpp_cvar <- function(x,type=1,var_names=NULL,save=FALSE,height=13,width=13,...)
{
    .plotbvar(x,type,var_names,save,height,width)
}

plot.Rcpp_dsge_gensys <- function(x,par_names=NULL,BinDenom=40,trace_plot=FALSE,save=FALSE,height=13,width=13,...)
{
    .plotedsge(x,par_names,BinDenom,trace_plot,save,height,width)
}

plot.Rcpp_dsge_uhlig <- function(x,par_names=NULL,BinDenom=40,trace_plot=FALSE,save=FALSE,height=13,width=13,...)
{
    .plotedsge(x,par_names,BinDenom,trace_plot,save,height,width)
}

plot.Rcpp_dsgevar_gensys <- function(x,par_names=NULL,BinDenom=40,MCMCplot=FALSE,save=FALSE,height=13,width=13,...)
{
    .plotdsgevar(x,par_names,BinDenom,MCMCplot,save,height,width)
}

plot.Rcpp_dsgevar_uhlig <- function(x,par_names=NULL,BinDenom=40,MCMCplot=FALSE,save=FALSE,height=13,width=13,...)
{
    .plotdsgevar(x,par_names,BinDenom,MCMCplot,save,height,width)
}

#

.plotbvar <- function(obj,type=1,var_names=NULL,save=FALSE,height=13,width=13)
{    
    obj_class <- class(obj)[1]

    constant <- obj$cons_term
    p <- obj$p
    c_int = obj$c_int

    if (obj_class == "Rcpp_bvars") {
        c_int = 0
    }
    
    K <- dim(obj$beta_draws)[1]
    M <- dim(obj$beta_draws)[2]
    n_draws <- dim(obj$beta_draws)[3]

    #
    
    BetaPerm <- aperm(obj$beta_draws,c(3,1,2))
    
    CoefLabels<-character(p*M)
    jj <- 1
    for(i in 1:p){
        for(k in 1:M){
            CoefLabels[jj]<-paste("L",i,".",k,sep="")
            jj = jj + 1
        }
    }

    if (class(var_names) != "character") {
        var_names <- character(length=M)
        for (i in 1:M) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }

    plotSigma <- FALSE
    if (obj_class == "Rcpp_bvars" || obj_class == "Rcpp_bvarw" || obj_class == "Rcpp_cvar") {
        plotSigma <- TRUE
    }

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    
    #

    CoefCount <- 1
    BinDenom <- 40
    ParamBin <- 1

    if (obj_class == "Rcpp_bvars") {
        q <- dim(obj$Psi_draws)[1]
        PsiPerm <- aperm(obj$Psi_draws,c(3,1,2))

        if(save==TRUE){
            if(class(dev.list()) != "NULL"){dev.off()}
            cairo_ps(filename="Psi.eps",height=(floor(height*q/M)),width=width)
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(q,M)))
        
        for(i in 1:q){
            for(j in 1:M){
                VarName <- var_names[j]
                #
                CFDF <- data.frame(PsiPerm[,i,j])
                colnames(CFDF) <- "CFDF"
                #
                if(j==1){
                    if(type==1){
                        ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                        print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Psi") + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
                    }else{
                        print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Psi") + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
                    }
                }else{
                    if(type==1){
                        ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                        print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
                    }else{
                        print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
                    }
                }
                Sys.sleep(0.3)
            }
            if(save==TRUE){dev.off()}
        }
    } else if (constant==TRUE) {
        
        if(save==TRUE){
            if(class(dev.list()) != "NULL"){dev.off()}
            cairo_ps(filename="Constant.eps",height=(floor(height/M)),width=width)
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(1,M)))
        
        for(j in 1:M){
            VarName <- var_names[j]
                
            CFDF <- data.frame(BetaPerm[,i,j])
            colnames(CFDF) <- "CFDF"
                
            if(j==1){
                if (type==1) {
                    ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                    print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Constant") + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
                } else {
                    print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab("Constant") + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
                }
            } else {
                if (type==1) {
                    ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                    print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkred",fill="black",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(1,j))
                } else {
                    print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkred",fill="red",alpha=0.25),vp = vplayout(1,j))
                }
            }
            Sys.sleep(0.3)
        }
        if(save==TRUE){dev.off()}
    }

    #

    for(i in 1:p){
        
        if(save==TRUE){
            if(class(dev.list()) != "NULL"){dev.off()}

            SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
            cairo_ps(filename=SaveLag,height=height,width=width)
        }
        
        #
        
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(M,M)))
        
        #
        
        for (j in 1:M) {
            for (l in 1:M) {
                VarName <- var_names[l]
                #
                CFDF <- data.frame(BetaPerm[,((i-1)*M+j+c_int),l])
                colnames(CFDF) <- "CFDF"
                #
                #j==1 is for the variable title;l==1 is the coefficient label on y-axis
                if (j==1) {
                    if (l==1) {
                        if (type==1) {
                            ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                            CoefCount <- CoefCount + 1
                        } else {
                            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                            CoefCount <- CoefCount + 1
                        }
                    } else {
                        if (type==1) {
                            ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=paste(VarName)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                        }else{
                            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarName)) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                        }
                    }
                } else {
                    if (l==1) {
                        if (type==1) {
                            ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                            CoefCount <- CoefCount + 1
                        } else {
                            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(paste(CoefLabels[CoefCount])) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                            CoefCount <- CoefCount + 1
                        }
                    } else {
                        if (type==1) {
                            ParamBin <- (max(CFDF) - min(CFDF))/BinDenom
                            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="darkblue",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                        } else {
                            print(ggplot(CFDF,aes(x=CFDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="blue",alpha=0.25),vp = vplayout(j,l))
                        }
                    }
                }
                Sys.sleep(0.3)
            }
        }
        if(save==TRUE){dev.off()}
    }
    
    #
    
    if (plotSigma==TRUE) {
        Sigmas <- obj$Sigma_draws
        SigmaPerm <- aperm(Sigmas,c(3,1,2))
        
        #

        if (save==TRUE) {
            if(class(dev.list()) != "NULL"){dev.off()}

            SaveSig <- paste("Sigma.eps",sep="")
            cairo_ps(filename=SaveSig,height=height,width=width)
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(M,M)))
        
        #

        for(j in 1:M){
            for(l in 1:M){
                VarNameY <- var_names[j]
                VarNameX <- var_names[l]
                #
                SDF <- data.frame(SigmaPerm[,j,l])
                colnames(SDF) <- "SDF"
                #
                #j==1 is for the variable title;l==1 is the coefficient label on y-axis
                if(j==1){
                    if(l==1){
                        if(type==1){
                            ParamBin <- (max(SDF) - min(SDF))/BinDenom
                            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + labs(title=paste(VarNameX)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                        }else{
                            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
                        }
                    }else{
                        if(type==1){
                            ParamBin <- (max(SDF) - min(SDF))/BinDenom
                            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + labs(title=paste(VarNameX)) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                        }else{
                            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + labs(title=paste(VarNameX)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
                        }
                    }
                }else{
                    if(l==1){
                        if(type==1){
                            ParamBin <- (max(SDF) - min(SDF))/BinDenom
                            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                        }else{
                            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(paste(VarNameY)) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
                        }
                    }else{
                        if(type==1){
                            ParamBin <- (max(SDF) - min(SDF))/BinDenom
                            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_histogram(colour="orchid4",fill="black",binwidth=ParamBin) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(j,l))
                        }else{
                            print(ggplot(SDF,aes(x=SDF)) + xlab(NULL) + ylab(NULL) + geom_density(colour="darkblue",fill="orchid4",alpha=0.25),vp = vplayout(j,l))
                        }
                    }
                }
                Sys.sleep(0.3)
            }
        }
        if(save==TRUE){dev.off()}
    }
}

.plotbvartvp <- function(obj,var_names=NULL,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13)
{
    obj_class <- class(obj)[1]

    constant <- obj$cons_term
    c_int = obj$c_int
    p <- obj$p
    K <- obj$K
    M <- obj$M
    n_draws <- dim(obj$alpha_draws)[3]
    kT <- dim(obj$alpha_draws)[2]

    #

    alpha_draws <- obj$alpha_draws

    alpha_draws <- aperm(alpha_draws,c(3,2,1)) # n_draws x (n - tau) x (K*M)
    alpha_draws_sorted <- apply(alpha_draws,c(2,3),sort)

    #

    CTPUpper <- round(percentiles[3]*n_draws)
    CTPMid <- round(percentiles[2]*n_draws)
    CTPLower <- round(percentiles[1]*n_draws)

    # 
    
    CoefLabels<-character(p*M)
    jj <- 1
    for(i in 1:p){
        for(k in 1:M){
            CoefLabels[jj]<-paste("L",i,".",k,sep="")
            jj = jj + 1
        }
    }

    if (class(var_names) != "character") {
        var_names <- character(length=M)
        for (i in 1:M) {  
            var_names[i] <- paste("VAR",i,sep="")
        }
    }

    #

    coefs_index <- matrix(1:(K*M),K,M,byrow=TRUE)

    for (i in 1:p) {
        coefs_index[(1+c_int+(i-1)*M):(c_int+i*M),] <- t(coefs_index[(1+c_int+(i-1)*M):(c_int+i*M),])
    }

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    
    if(class(dev.list()) != "NULL"){dev.off()}

    Time <- CTPL <- CTPM <- CTPU <- NULL
    CoefCount <- 1
    
    #

    if (constant==TRUE) {
        for (i in 1:1) {
            if(save==TRUE){
                cairo_ps(filename="Constant.eps",height=(floor(height/M)),width=width)
            }

            pushViewport(viewport(layout=grid.layout(1,M)))
            
            for (j in 1:M) {
                VarName <- var_names[j]

                CFDF <- data.frame(alpha_draws_sorted[CTPLower,,j],alpha_draws_sorted[CTPMid,,j],alpha_draws_sorted[CTPUpper,,j],tau:(kT+tau-1))
                colnames(CFDF) <- c("CTPL","CTPM","CTPU","Time")

                if (j==1) {
                    print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab("Constant") + labs(title=paste(VarName)) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(1,j))
                } else {
                    print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(NULL) + labs(title=paste(VarName)) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(1,j))
                }

                Sys.sleep(0.3)
            }
            if(save==TRUE){dev.off()}
        }
    }

    #

    for (i in 1:p) {
        
        if(save==TRUE){
            SaveLag <- paste("CoefLag",as.character(i),".eps",sep="")
            cairo_ps(filename=SaveLag,height=height,width=width)
        }
        
        grid.newpage()
        pushViewport(viewport(layout=grid.layout(M,M)))
        
        #

        for (j in 1:M) {
            for (l in 1:M) {
                VarName <- var_names[l]
                coef_ind <- coefs_index[c_int+(i-1)*M+j,l]

                CFDF <- data.frame(alpha_draws_sorted[CTPLower,,coef_ind],alpha_draws_sorted[CTPMid,,coef_ind],alpha_draws_sorted[CTPUpper,,coef_ind],tau:(kT+tau-1))
                colnames(CFDF) <- c("CTPL","CTPM","CTPU","Time")

                # j==1 is for the variable title; l==1 is the coefficient label on y-axis

                if (j==1) {
                    if (l==1) {
                        print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(paste(CoefLabels[CoefCount])) + labs(title=paste(VarName)) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1)+ geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(j,l))
                        CoefCount <- CoefCount + 1
                    } else {
                        print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(NULL) + labs(title=paste(VarName)) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(j,l))
                    }
                } else {
                    if (l==1) {
                        print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(paste(CoefLabels[CoefCount])) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(j,l))
                        CoefCount <- CoefCount + 1
                    } else {
                        print(ggplot(data=(CFDF),aes(x=Time)) + xlab("Time") + ylab(NULL) + geom_ribbon(aes(ymin=CTPL,ymax=CTPU),color="lightslateblue",lty=1,fill="lightslateblue",alpha=0.2,size=0.1) + geom_line(aes(y=CTPM),color="darkgreen",size=2),vp = vplayout(j,l))
                    }
                }
                Sys.sleep(0.3)
            }
        }

        #

        if(save==TRUE){dev.off()}
    }
    #
}

.plotedsge <- function(obj,par_names=NULL,BinDenom=40,trace_plot=FALSE,save=FALSE,height=13,width=13)
{
    n_draws <- dim(obj$dsge_draws)[1]
    n_param <- dim(obj$dsge_draws)[2]

    if (n_draws==0) {
        stop("no MCMC draws detected.\n",call.=FALSE)
    }

    dsge_draws <- obj$dsge_draws # make a copy

    #

    vplayout <- function(x,y){viewport(layout.pos.row=x, layout.pos.col=y)}
    
    #
    
    if (is.null(par_names)==TRUE) {
        par_names <- character(length=n_param)
        for(i in 1:n_param){  
            par_names[i] <- paste("parameter",i,sep="")
        }
    }
    
    #

    MR <- 0; MC <- 0
    plot_pages <- 1
    
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

    param_count <- 1

    for (j in 1:plot_pages) {
        
        if (save==TRUE) {
            if(class(dev.list()) != "NULL"){dev.off()}
            
            if (plot_pages==1) {
                cairo_ps(filename="DSGEParameters.eps",height=height,width=width)
            } else {
                SaveParam <- paste("DSGEParameters_",j,".eps",sep="")
                cairo_ps(filename=SaveParam,height=height,width=width)
            }
        }

        grid.newpage()
        pushViewport(viewport(layout=grid.layout(MR,MC)))

        #

        DSGEPar <- Keep <- NULL # CRAN check workaround
        
        for (i in 1:MR) {
            for (k in 1:MC) {
                #
                if (param_count <= n_param) {

                    ParamDF <- data.frame(dsge_draws[,param_count])
                    colnames(ParamDF) <- c("DSGEPar")
                    
                    #

                    ParamBin <- (max(ParamDF) - min(ParamDF))/BinDenom
                    print(ggplot(data=(ParamDF),aes(DSGEPar)) + xlab("") + ylab("") + geom_histogram(colour="darkblue",binwidth=ParamBin) + labs(title=par_names[param_count]) + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                    
                    #

                    param_count <- param_count + 1
                    Sys.sleep(0.3)

                } else {param_count <- param_count + 1}
            }
        }

        if(save==TRUE){dev.off()}
    }

    #
    # Trace Plot

    if (trace_plot==TRUE) {

        param_count <- 1

        for (j in 1:plot_pages) {
            #
            if (save==TRUE) {
                if(class(dev.list()) != "NULL"){dev.off()}

                if (plot_pages==1) {
                    cairo_ps(filename="MCMCplot.eps",height=height,width=width)
                } else {
                    SaveParam <- paste("MCMCplot_",j,".eps",sep="")
                    cairo_ps(filename=SaveParam,height=height,width=width)
                }
            }
            
            grid.newpage()
            pushViewport(viewport(layout=grid.layout(MR,MC)))
            
            #
            
            for(i in 1:MR){
                for(k in 1:MC){
                    #
                    if(param_count <= n_param){
                        ParamDF <- data.frame(dsge_draws[,param_count],1:n_draws)
                        colnames(ParamDF) <- c("DSGEPar","Keep")
                        
                        #
                        
                        print(ggplot(data=(ParamDF),aes(x=Keep)) + xlab("Keep Run") + ylab(paste(par_names[param_count])) + geom_line(aes(y=DSGEPar),color="black") + theme(panel.background = element_rect(fill='white', colour='grey5')) + theme(panel.grid.major = element_line(colour = 'grey89')),vp = vplayout(i,k))
                        
                        #
                        
                        param_count <- param_count + 1
                        Sys.sleep(0.3)
                        
                    }else{param_count <- param_count + 1}
                }
            }
            
            if(save==TRUE){dev.off()}
        }
    }
    #
}

.plotdsgevar <- function(obj,par_names=NULL,BinDenom=40,trace_plot=FALSE,save=FALSE,height=13,width=13)
{    
    dsge_obj <- obj$dsge
    .plotedsge(dsge_obj,par_names,BinDenom,trace_plot,save,height,width)
}
