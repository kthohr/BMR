# 01/23/15
prior <- function(priorform,priorpars,parname=NULL,moments=TRUE,NR=NULL,NC=NULL){
  .priorfunction(priorform,priorpars,parname,moments,NR,NC)
}

.priorfunction <- function(priorform,priorpars,parname=NULL,moments=TRUE,NR=NULL,NC=NULL){
  #
  x <- 0; PriorDF <- 0;
  xmean <- 0; xmode <- 0; xvar <- 0
  #
  if(priorform=="Beta"){
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
    if(moments==TRUE){
      cat('Mean: ', xmean, '. \n', sep="")
      cat('Mode: ', xmode, '. \n', sep="")
      cat('Variance: ', xvar, '. \n', sep="")
    }
    #
    PriorDF <- data.frame(hx,x)
    colnames(PriorDF) <- c("Density","Domain")
    #
  }else if(priorform=="Uniform"){
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
  }else if(priorform=="IGamma"){
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
  }else if(priorform=="Gamma"){
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
  }else if(priorform=="Normal"){
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
  }
  #
  #
  #
  if(class(NR)=="NULL" | class(NC)=="NULL"){
    if(class(parname)=="character"){
      if(priorform=="Uniform"){
        print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_line(aes(y=Density),color="red",size=1) + labs(title=parname))
      }else{
        print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_vline(xintercept=xmode,linetype = "longdash") + geom_line(aes(y=Density),color="red",size=1) + labs(title=parname))
      }
    }else{
      if(priorform=="Uniform"){
        print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_line(aes(y=Density),color="red",size=1))
      }else{
        print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_vline(xintercept=xmode,linetype = "longdash") + geom_line(aes(y=Density),color="red",size=1))
      }
    }
  }else{
    if(class(parname)=="character"){
      if(priorform=="Uniform"){
        print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_line(aes(y=Density),color="red",size=1) + labs(title=parname),vp = vplayout(NR,NC))
      }else{
        print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_vline(xintercept=xmode,linetype = "longdash") + geom_line(aes(y=Density),color="red",size=1) + labs(title=parname),vp = vplayout(NR,NC))
      }
    }else{
      if(priorform=="Uniform"){
        print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_line(aes(y=Density),color="red",size=1),vp = vplayout(NR,NC))
      }else{
        print(ggplot(data=(PriorDF),aes(x=Domain)) + xlab("") + geom_ribbon(aes(ymin=0,ymax=Density),color="blue",lty=1,fill="blue",alpha=0.2,size=0.1) + geom_vline(xintercept=xmean) + geom_vline(xintercept=xmode,linetype = "longdash") + geom_line(aes(y=Density),color="red",size=1),vp = vplayout(NR,NC))
      }
    }
  }
}