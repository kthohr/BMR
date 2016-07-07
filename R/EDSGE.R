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

EDSGE.default <- function(dsgedata,chains=1,cores=1,
                          ObserveMat,initialvals,partomats,
                          priorform,priorpars,parbounds,parnames=NULL,
                          optimMethod="Nelder-Mead",
                          optimLower=NULL,optimUpper=NULL,
                          optimControl=list(),
                          DSGEIRFs=TRUE,irf.periods=20,
                          scalepar=1,keep=50000,burnin=10000,
                          tables=TRUE){
    #
    message('Trying to solve the model with your initial values... ',appendLF=FALSE)
    dsgemats1t <- partomats(initialvals)
    dsgesolved1t <- SDSGE(dsgemats1t)
    StateMats1t <- statespace(dsgesolved1t)
    soltype <- NULL
    if(dsgesolved1t$sol_type==1){
        soltype <- "gensys."
    }else if(dsgesolved1t$sol_type==2){
        soltype <- "Uhlig\'s method."
    }
    message('Done, using ', soltype)
    #
    priorformRet <- priorform
    prelimwork <- .edsgePrelimWork(dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
    priorform <- prelimwork$priorform; parbounds <- prelimwork$parbounds
    #
    #
    #
    Mode <- .DSGEModeEstimate(dsgedata,ObserveMat,initialvals,partomats,priorform,priorpars,parbounds,parnames,optimMethod,optimLower,optimUpper,optimControl,tables)
    #
    dsgemode <- Mode$dsgemode; parMode <- Mode$parMode; parModeSEs <- Mode$parModeSEs
    #
    dsgeret <- 0
    if(keep==0){
        parRet <- matrix(0,nrow=0,ncol=length(parMode))
        if(class(parnames)=="character"){
            colnames(parRet) <- parnames
        }
        dsgeret <- list(Parameters=parRet,parMode=parMode,IRFs=NULL,ModeHessian=dsgemode$hessian,logMargLikelihood=Mode$logMargLikelihood,scalepar=scalepar,AcceptanceRate=NULL,RootRConvStats=NULL,ObserveMat=ObserveMat,data=dsgedata,partomats=partomats,priorform=priorformRet,priorpars=priorpars,parbounds=parbounds)
        class(dsgeret) <- "EDSGE"
        #
        return(dsgeret)
    }
    #
    #
    #
    message(' ', sep="")
    message('Beginning MCMC run, ', date(),'.', sep="")
    MCMCRes <- 0
    if(chains==1){
        MCMCRes <- .DSGEMCMC(dsgemode,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
    }else if(chains > 1){
        MCMCRes <- .DSGEMCMCMulti(dsgemode,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds,chains,cores)
    }
    message('MCMC run finished, ', date(),'.', sep="")
    #
    if(class(parnames)=="character"){
        colnames(MCMCRes$parameters) <- parnames
    }
    #
    PostMCMCInfo <- .DSGEMCMCPrint(MCMCRes,chains,parMode,parModeSEs,parnames,tables)
    #
    #
    #
    IRFs <- NULL
    if(DSGEIRFs == TRUE){
        message(' ', sep="")
        message('Computing IRFs now... ',appendLF=FALSE)
        IRFs <- array(0,dim=c(irf.periods,ncol(StateMats1t$F),ncol(StateMats1t$G),keep))
        for(i in 1:keep){
            dsgemats <- partomats(MCMCRes$parameters[i,])
            dsgesolved <- SDSGE(dsgemats)
            iIRF <- IRF(dsgesolved,sqrt(diag(dsgemats$shocks)),irf.periods,varnames=NULL,plot=FALSE,save=FALSE)
            IRFs[,,,i] <- iIRF$IRFs
        }
        message('Done.')
    }
    #
    dsgeret <- list(Parameters=MCMCRes$parameters,parMode=parMode,ModeHessian=dsgemode$hessian,logMargLikelihood=Mode$logMargLikelihood,IRFs=IRFs,scalepar=scalepar,AcceptanceRate=MCMCRes$acceptRate,RootRConvStats=PostMCMCInfo$Diagnostics,ObserveMat=ObserveMat,data=dsgedata,partomats=partomats,priorform=priorformRet,priorpars=priorpars,parbounds=parbounds)
    class(dsgeret) <- "EDSGE"
    #
    return(dsgeret)
}

.edsgePrelimWork <- function(dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds){
    #
    # Change from character to numeric
    priorformNum <- numeric(length(priorform))
    # Normal = 1, Gamma = 2, IGamma = 3, Beta = 4, Uniform = 5
    for(i in 1:length(priorform)){
        if(priorform[i]=="Normal"){
            priorformNum[i] <- 1
        }else if(priorform[i]=="Gamma"){
            priorformNum[i] <- 2
        }else if(priorform[i]=="IGamma"){
            priorformNum[i] <- 3
        }else if(priorform[i]=="Beta"){
            priorformNum[i] <- 4
        }else if(priorform[i]=="Uniform"){
            priorformNum[i] <- 5
        }else{
            stop("Parameter ", i ," does not have a valid prior form.\n",call.=FALSE)
        }
    }
    #
    # Check if parbounds is set correctly for uniform priors
    for(i in 1:length(priorform)){
        if(priorform[i]=="Uniform"){
            parbounds[i,] <- priorpars[i,]
        }
    }
    #
    return=list(priorform=priorformNum,parbounds=parbounds)
}

.DSGEParTransform <- function(parameters,priorform,parbounds,trans=1){
    #
    if(trans == 1){
        nparam <- as.numeric(length(parameters))
        parametersTrans <- numeric(nparam)
        #
        for(i in 1:nparam){
            if(priorform[i] == 2 || priorform[i] == 3){
                if(is.na(parbounds[i,1]) != TRUE){
                    parametersTrans[i] <- log(parameters[i] - parbounds[i,1])
                }else{
                    parametersTrans[i] <- log(parameters[i])
                }
            }else if(priorform[i] == 4 || priorform[i] == 5){
                parametersTrans[i] <- log((parameters[i] - parbounds[i,1])/(parbounds[i,2] - parameters[i]))
            }else{#Normal
                parametersTrans[i] <- parameters[i]
            }
        }
        #
        return(parametersTrans)
    }else{
        nparam <- as.numeric(length(parameters))
        TransBack <- numeric(nparam)
        #
        for(i in 1:nparam){
            if(priorform[i] == 2 || priorform[i] == 3){
                if(is.na(parbounds[i,1]) != TRUE){
                    TransBack[i] <- exp(parameters[i]) + parbounds[i,1]
                }else{
                    TransBack[i] <- exp(parameters[i])
                }
            }else if(priorform[i] == 4 || priorform[i] == 5){
                TransBack[i] <- (parbounds[i,1] + parbounds[i,2]*exp(parameters[i]))/(1 + exp(parameters[i]))
            }else{#Normal
                TransBack[i] <- parameters[i]
            }
        }
        #
        return(TransBack)
    }
}

.DSGEPriors <- function(parameters,parametersTrans,priorform,priorpars,parbounds,dsgelike){
    #
    nparam <- as.numeric(length(parameters))
    dsgeposterior <- - dsgelike
    #
    for(i in 1:nparam){
        if(priorform[i] == 1){
            dsgeposterior <- dsgeposterior + dnorm(parameters[i],mean = priorpars[i,1],sd = sqrt(priorpars[i,2]), log = TRUE)
        }else if(priorform[i] == 2){
            dsgeposterior <- dsgeposterior + dgamma(parameters[i],shape=priorpars[i,1],scale=priorpars[i,2],log=TRUE) + parametersTrans[i]
        }else if(priorform[i] == 3){
            dsgeposterior <- dsgeposterior + log( ( (priorpars[i,2]^priorpars[i,1])/gamma(priorpars[i,1]) )*( (parameters[i])^(-priorpars[i,1]-1) )*( exp(-priorpars[i,2]/parameters[i]) ) ) + parametersTrans[i]
        }else if(priorform[i] == 4){
            if(is.na(parbounds[i,1]) || is.na(parbounds[i,2])){
                z = (parameters[i]-parbounds[i,1])/(parbounds[i,2]-parbounds[i,1])
                dsgeposterior <- dsgeposterior + dbeta(z,priorpars[i,1],priorpars[i,2],log=TRUE)
                dsgeposterior <- dsgeposterior + parametersTrans[i] - 2*log(1+exp(parametersTrans[i]))
            }else{
                z = (parameters[i]-parbounds[i,1])/(parbounds[i,2]-parbounds[i,1])
                dsgeposterior <- dsgeposterior + dbeta(z,priorpars[i,1],priorpars[i,2],log=TRUE)
                dsgeposterior <- dsgeposterior + log(parbounds[i,2] - parbounds[i,1]) + parametersTrans[i] - 2*log(1+exp(parametersTrans[i]))
            }
        }else if(priorform[i]==5){
            dsgeposterior <- dsgeposterior + dunif(parameters[i],priorpars[i,1],priorpars[i,2],log=TRUE)
            dsgeposterior <- dsgeposterior + log(parbounds[i,2] - parbounds[i,1]) + parametersTrans[i] - 2*log(1+exp(parametersTrans[i]))
        }else{
            stop("Parameter ", i ," does not have a valid prior form.\n",call.=FALSE)
        }
    }
    # Return negative of the posterior:
    dsgeposterior <- - dsgeposterior
    #
    return(as.numeric(dsgeposterior))
}

.dsgeposteriorfn <- function(parametersTrans,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds){
    #
    parameters <- .DSGEParTransform(parametersTrans,priorform,parbounds,2)
    #
    dsgemats <- partomats(parameters)
    dsgesolved <- SDSGE(dsgemats)
    #
    StateMats <- statespace(dsgesolved)
    #
    #dsgelike <- .Call("DSGEKalman", dsgedata,ObserveMat,dsgemats$ObsCons,StateMats$F,StateMats$G,dsgemats$shocks,dsgemats$MeasErrs,200, PACKAGE = "BMR")
    dsgelike <- .Call("DSGECR", dsgedata,ObserveMat,dsgemats$ObsCons,StateMats$F,StateMats$G,dsgemats$shocks,dsgemats$MeasErrs,200, PACKAGE = "BMR")
    #
    dsgeposterior <- .DSGEPriors(parameters,parametersTrans,priorform,priorpars,parbounds,dsgelike$dsgelike)
    #
    return(dsgeposterior)
}

.LaplaceMargLikelihood <- function(obj){
    PosteriorVal <- -obj$value
    MargLike <- PosteriorVal + (length(obj$par)*log(2*pi) + log(1/det(obj$hessian)))/2
    return(MargLike)
}

.DSGEModeEstimate <- function(dsgedata,ObserveMat,initialvals,partomats,
                              priorform,priorpars,parbounds,parnames,
                              optimMethod,optimLower,optimUpper,optimControl,
                              tables){
    #
    parametersTrans <- .DSGEParTransform(initialvals,priorform,parbounds,1)
    #
    OptimMethods <- optimMethod
    #
    dsgemode <- NULL
    prevlogpost <- 0
    message(' ', sep="")
    message('Beginning optimization, ', date(),'.', sep="")
    for(jj in 1:length(OptimMethods)){
        #
        optimMethod <- OptimMethods[jj]
        #
        if(jj==1){
            message('Using Optimization Method: ',optimMethod,'. ', sep="")
        }else{
            message('Using Optimization Method: ',optimMethod,'. ', sep="",appendLF = FALSE)
        }
        #
        if(optimMethod=="Nelder-Mead"){
            dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="Nelder-Mead",control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
        }else if(optimMethod=="BFGS"){
            dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="BFGS",control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
        }else if(optimMethod=="CG"){
            dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="CG",control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
        }else if(optimMethod=="L-BFGS-B"){
            dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="L-BFGS-B",lower=optimLower,upper=optimUpper,control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
        }else if(optimMethod=="SANN"){
            dsgemode <- optim(parametersTrans,fn=.dsgeposteriorfn,method="SANN",control=optimControl,dsgedata=dsgedata,ObserveMat=ObserveMat,partomats=partomats,priorform=priorform,priorpars=priorpars,parbounds=parbounds,hessian=TRUE)
        }else{
            stop("You have entered an unrecognized optimization method.\n",call.=FALSE)
        }
        #
        if(jj>1){
            message('Change in the log posterior: ',round(-dsgemode$value - prevlogpost,5),'. ', sep="")
        }
        #
        parametersTrans <- dsgemode$par
        prevlogpost <- -dsgemode$value
    }
    #
    ConvCode <- dsgemode$convergence
    ConvReport <- 0
    if(ConvCode==0){
        ConvReport <- "successful completion"
    }else if(ConvCode==1){
        ConvReport <- "maximum number of iterations reached"
    }else if(ConvCode==10){
        ConvReport <- "degeneracy of the Nelder-Mead simplex"
    }else if(ConvCode==51){
        ConvReport <- "warning from L-BFGS-B"
    }else if(ConvCode==52){
        ConvReport <- "error from L-BFGS-B"
    }else{
        ConvReport <- "unknown"
    }
    #
    message('Optimization over, ', date(),'. ', sep="")
    message(' ', sep="")
    message('Optimizer Convergence Code: ',dsgemode$convergence,'; ',ConvReport,'. ', sep="")
    if(tables==TRUE){
        message(' ', sep="")
        cat('Optimizer Iterations: \n', sep="")
        print(dsgemode$counts)
    }
    #
    parMode <- .DSGEParTransform(dsgemode$par,priorform,parbounds,2)
    #
    logMargLikelihood <- .LaplaceMargLikelihood(dsgemode)
    #
    message(' ', sep="")
    message('Log Marginal Likelihood: ',round(logMargLikelihood,5),'.', sep="")
    #
    parModeHessian <- solve(dsgemode$hessian)
    parModeHessian <- diag(parModeHessian)
    parModeHessian <- sqrt(parModeHessian)
    #
    parModeSEs <- .DSGEParTransform(dsgemode$par,priorform,parbounds,2) - .DSGEParTransform(dsgemode$par - parModeHessian,priorform,parbounds,2)
    #
    parMode <- matrix(parMode,nrow=1)
    parModeSEs <- matrix(parModeSEs,nrow=1)
    parModeHessian <- matrix(parModeHessian,nrow=1)
    #
    ModeTable <- matrix(NA,nrow=length(dsgemode$par),ncol=2)
    ModeTable[,1] <- parMode
    ModeTable[,2] <- parModeSEs
    #
    colnames(ModeTable) <- c("Estimate","SE")
    if(class(parnames)=="character"){
        rownames(ModeTable) <- parnames
    }
    if(tables==TRUE){
        message(' ', sep="")
        message('Parameter Estimates and Standard Errors (SE) at the Posterior Mode: ', sep="")
        message(' ', sep="")
        print(ModeTable)
    }
    #
    rownames(parMode) <- "Parameter:"
    rownames(parModeSEs) <- "Parameter:"
    rownames(parModeHessian) <- "Parameter:"
    if(class(parnames)=="character"){
        colnames(parMode) <- parnames
        colnames(parModeHessian) <- parnames
        colnames(parModeSEs) <- parnames
    }
    #
    #
    return=list(dsgemode=dsgemode,parMode=parMode,parModeSEs=parModeSEs,logMargLikelihood=logMargLikelihood)
}

.DSGEMCMC <- function(dsgeopt,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds,parallel=FALSE){
    #
    Draws <- matrix(NA,nrow=(keep+1),ncol=length(dsgeopt$par))
    #
    InitialDraw <- dsgeopt$par
    if(parallel==TRUE){
        InitialDraw <- dsgeopt$par + runif(1,-1,1)*c(sqrt(diag(solve(dsgeopt$hessian))))
    }
    #
    PrevLP <- (-1)*dsgeopt$value
    if(parallel==TRUE){
        PrevLP <- (-1)*.dsgeposteriorfn(c(InitialDraw),dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
    }
    #
    PickMeInstead <- matrix(c(InitialDraw))
    #
    CovM <- scalepar*(solve(dsgeopt$hessian))
    CovMChol <- t(chol(CovM))
    #
    Acceptances <- 0
    #
    for (i in 1:burnin){
        #
        proposal <- PickMeInstead + CovMChol%*%matrix(rnorm(length(dsgeopt$par)))
        #
        PropLP <- (-1)*.dsgeposteriorfn(c(proposal),dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
        if(is.nan(PropLP)){
            PropLP <- -1000000
        }
        #
        if(runif(1) < exp(PropLP-PrevLP)){
            PickMeInstead <- proposal
            PrevLP <- PropLP
        }
    }
    #
    Draws[1,] <- t(PickMeInstead)
    #
    for (i in 1:keep){
        #
        proposal <- PickMeInstead + CovMChol%*%matrix(rnorm(length(dsgeopt$par)))
        #
        PropLP <- (-1)*.dsgeposteriorfn(c(proposal),dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds)
        if(is.nan(PropLP)){
            PropLP <- -1000000
        }
        #
        if(runif(1) < exp(PropLP-PrevLP)){
            Draws[i+1,] <- t(proposal)
            Acceptances <- Acceptances + 1
            #
            PickMeInstead <- proposal
            PrevLP <- PropLP
        }else{
            Draws[i+1,] <- Draws[i,]
        }
    }
    #
    Draws <- Draws[-1,]
    accept <- Acceptances/keep
    #
    for(i in 1:keep){
        Draws[i,] <- .DSGEParTransform(Draws[i,],priorform,parbounds,2)
    }
    #
    if(parallel==FALSE){
        return=list(parameters=Draws,acceptRate=accept)
    }else{
        return(list(Draws,accept))
    }
}

.DSGEMCMCMulti <- function(dsgeopt,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds,chains,cores){
    #
    cl <- makeCluster(cores)
    registerDoParallel(cl)
    parallelsol <- 0
    parallelsol <- foreach(jj=1:chains, .packages=c("BMR")) %dopar% {
        .DSGEMCMC(dsgeopt,scalepar,keep,burnin,dsgedata,ObserveMat,partomats,priorform,priorpars,parbounds,TRUE)
    }
    #
    stopCluster(cl)
    #
    Draws <- matrix(NA,nrow=keep*chains,ncol=length(dsgeopt$par))
    accept <- numeric(chains)
    #
    for(j in 1:chains){
        Draws[((j-1)*keep+1):(j*keep),] <- parallelsol[[j]][[1]]
        accept[j] <- parallelsol[[j]][[2]]
    }
    #
    return=list(parameters=Draws,acceptRate=accept)
}

.MCMCDiagnostics <- function(parameters,chains){
    #
    keep <- nrow(parameters)/chains
    parArray <- array(0,dim=c(nrow(parameters),ncol(parameters),chains))
    for(i in 1:chains){
        parArray[,,i] <- parameters[((i-1)*keep+1):(i*keep),]
    }
    #
    meansAC <- numeric(dim(parArray)[2])
    for(i in 1:dim(parArray)[3]){
        for(j in 1:dim(parArray)[2]){
            meansAC[j] <- meansAC[j] + mean(parArray[,j,i])
        }
    }
    #
    meansAC <- meansAC/dim(parArray)[3]
    #
    B <- numeric(dim(parArray)[2])
    for(i in 1:dim(parArray)[3]){
        for(j in 1:dim(parArray)[2]){
            B[j] <- B[j] + (mean(parArray[,j,i]) - meansAC[j])^2
        }
    }
    #
    B <- dim(parArray)[1]/(dim(parArray)[3] - 1)*B
    #
    S <- numeric(dim(parArray)[2])
    for(i in 1:dim(parArray)[3]){
        for(j in 1:dim(parArray)[2]){
            S[j] <- S[j] + var(parArray[,j,i])
        }
    }
    W <- S/dim(parArray)[3]
    #
    #
    #
    rootR <- sqrt( ( ((dim(parArray)[1]-1)/dim(parArray)[1])*W + (1/dim(parArray)[1])*B ) / W)
    #
    return(rootR)
}

.DSGEMCMCPrint <- function(MCMCRes,chains,parMode,parModeSEs,parnames,tables){
    #
    Diagnostics <- NULL
    if(chains==1){
        message('Acceptance Rate: ', MCMCRes$acceptRate,'. ', sep="")
    }else{
        message('Acceptance Rate: ', sep="", appendLF=FALSE)
        for(kk in 1:(chains-1)){
            message('Chain ',kk,': ', MCMCRes$acceptRate[kk],'; ', sep="", appendLF=FALSE)
        }
        message('Chain ',chains,': ', MCMCRes$acceptRate[chains],'. ', sep="")
        #
        # Chain convergence statistics:
        #
        Diagnostics <- matrix(.MCMCDiagnostics(MCMCRes$parameters,chains),nrow=1)
        #
        rownames(Diagnostics) <- "Stat:"
        if(class(parnames)=="character"){
            colnames(Diagnostics) <- parnames
        }
        #
        if(tables==TRUE){
            message(' ', sep="")
            message('Root-R Chain-Convergence Statistics: ', sep="")
            print(Diagnostics)
        }
        message(' ', sep="")
    }
    #
    PostTable <- matrix(NA,nrow=length(parMode),ncol=4)
    PostTable[,1] <- parMode
    PostTable[,2] <- parModeSEs
    PostTable[,3] <- apply(MCMCRes$parameters,2,mean)
    PostTable[,4] <- apply(MCMCRes$parameters,2,sd)
    #
    colnames(PostTable) <- c("Posterior.Mode","SE.Mode","Posterior.Mean","SE.Posterior")
    if(class(parnames)=="character"){
        rownames(PostTable) <- parnames
    }
    if(tables==TRUE){
        message(' ', sep="")
        message('Parameter Estimates and Standard Errors: ', sep="")
        message(' ', sep="")
        print(PostTable)
    }
    #
    return=list(Diagnostics=Diagnostics)
}