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

DSGESim.SDSGE <- function(obj,shocks.cov,sim.periods,burnin=NULL,seedval=1122,hpfiltered=FALSE,lambda=1600,...){
  results <- .dsgesimulation(obj,shocks.cov,sim.periods,burnin,seedval,hpfiltered,lambda)
  return(results)
}

DSGESim.gensys <- function(obj,shocks.cov,sim.periods,burnin=NULL,seedval=1122,hpfiltered=FALSE,lambda=1600,...){
  results <- .dsgesimulation(obj,shocks.cov,sim.periods,burnin,seedval,hpfiltered,lambda)
  return(results)
}

DSGESim.uhlig <- function(obj,shocks.cov,sim.periods,burnin=NULL,seedval=1122,hpfiltered=FALSE,lambda=1600,...){
  results <- .dsgesimulation(obj,shocks.cov,sim.periods,burnin,seedval,hpfiltered,lambda)
  return(results)
}

.dsgesimulation <- function(obj,shocks.cov,sim.periods,burnin=NULL,seedval,hpfiltered=FALSE,lambda=1600){
  #
  StateMats <- statespace(obj)
  F <- StateMats$F; G <- StateMats$G
  #
  nShocks <- ncol(G)
  #
  if(class(burnin) != "numeric"){
    burnin <- ceiling(0.5*sim.periods)
  }
  #
  DSim <- matrix(0,nrow=(sim.periods+burnin),ncol=ncol(F))
  #
  set.seed(seedval)
  #
  A <- t(chol(shocks.cov))
  krun <- (sim.periods + burnin)*nShocks
  samp <- matrix(rnorm(krun),ncol=nShocks)
  shocks <- t(samp%*%t(A))
  #
  DSim[1,] <- t(G%*%shocks[,1])
  for(i in 2:(burnin+sim.periods)){
    DSim[i,] <- DSim[i-1,]%*%t(F) + t(G%*%shocks[,i])
  }
  #
  DSim <- DSim[(burnin+1):(burnin+sim.periods),]
  #
  if(hpfiltered==TRUE){
    hpfilterq <- function(x,lambda=1600){
      eye <- diag(length(x))
      result <- solve(eye+lambda*crossprod(diff(eye,lag=1,d=2)),x)
      return(result)
    }
    for(j in 1:ncol(DSim)){
      DSim[,j]<-hpfilterq(DSim[,j])
    }
  }
  #
  return(DSim)
}