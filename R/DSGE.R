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

SDSGE.default <-function(mats,type=NULL){
  #
  sdsgerest <- 0
  #
  if(is.null(type)){ # type=NULL
    if(is.null(mats$Gamma0) && is.null(mats$A)){
      stop("Unrecognized solution method. Check that your matrices are correctly labelled.\n",call.=FALSE)
    }
    # gensys
    if((is.null(mats$Gamma0)==FALSE) && (is.null(mats$A)==TRUE)){
      sdsgeret <- gensys(mats$Gamma0, mats$Gamma1, mats$C, mats$Psi, mats$Pi)
      sdsgeret$sol_type <- 1
    }
    # uhlig
    if((is.null(mats$Gamma0)==TRUE) && (is.null(mats$A)==FALSE)){
      sdsgeret <- uhlig(mats$A,mats$B,mats$C,mats$D,mats$F,mats$G,mats$H,mats$J,mats$K,mats$L,mats$M,mats$N,mats$whichEig)
      sdsgeret$sol_type <- 2
    }
  }else if(type==1){ #gensys
    if(is.null(mats$Gamma0)){
      stop("Gensys mats not labelled correctly.\n",call.=FALSE)
    }else{
      sdsgeret <- gensys(mats$Gamma0, mats$Gamma1, mats$C, mats$Psi, mats$Pi)
      sdsgeret$sol_type <- 1
    }
  }else if(type==2){ #uhlig
    if(is.null(mats$A)){
      stop("Uhlig mats not labelled correctly.\n",call.=FALSE)
    }else{
      sdsgeret <- uhlig(mats$A,mats$B,mats$C,mats$D,mats$F,mats$G,mats$H,mats$J,mats$K,mats$L,mats$M,mats$N,mats$whichEig)
      sdsgeret$sol_type <- 2
    }
  }else{
    stop("Unrecognized solution type. type should be 1 (gensys) or 2 (uhlig).\n",call.=FALSE)
  }
  #
  class(sdsgeret) <- "SDSGE"
  return(sdsgeret)
}

gensys.default <- function(Gamma0,Gamma1,C,Psi,Pi){
  #
  SimsSol <- .Call("gensysCpp", Gamma0,Gamma1,C,Psi,Pi, PACKAGE = "BMR")
  #
  gensysret <- list(G1=SimsSol$G1,Cons=SimsSol$Cons,impact=SimsSol$impact,eu=SimsSol$eu,
                    Psi=Psi,Pi=Pi)
  #
  class(gensysret) <- "gensys"
  return(gensysret)
}

uhlig.default <-function(A,B,C,D,F,G,H,J,K,L,M,N,whichEig=NULL){
  #
  if(is.null(whichEig)){
    whichEig <- numeric(0)
  }else{
    whichEig <- c(whichEig)
  }
  #
  UhligSol <- .Call("UhligCpp", A,B,C,D,F,G,H,J,K,L,M,N,whichEig, PACKAGE = "BMR")
  #
  uhligret <- list(N=N,P=UhligSol$P,Q=UhligSol$Q,R=UhligSol$R,S=UhligSol$S,EigenValues=UhligSol$EigValueSorted,EigenVectors=UhligSol$EigVecSorted)
  class(uhligret) <- "uhlig"
  return(uhligret)
}

statespace.gensys <- function(obj){
  #
  Psi <- obj$Psi; Pi <- obj$Pi
  G1 <- obj$G1; impact <- obj$impact; Cons <- obj$Cons
  #
  n <- nrow(G1)
  #
  # we first need to shift some elements around:
  # drop expectational variables
  #
  nshocks <- ncol(Psi)
  #
  expec_ind <- matrix((1:n),nrow=1)%*%(Pi != 0)
  expec_ind <- as.numeric(unique(expec_ind))
  #
  indt <- (1:n)[-c(expec_ind)]
  #
  G1 <- G1[indt,indt] 
  if(nshocks==1){
    impact <- matrix(impact[indt])
  }else{
    impact <- impact[indt,]
  }
  Cons <- matrix(Cons[indt])
  #
  # Now form F and G
  #
  F <- G1
  G <- impact
  #
  return=list(F=F,G=G)
}

statespace.uhlig <- function(obj){
  #
  N <- obj$N
  P <- obj$P; Q <- obj$Q; R <- obj$R; S <- obj$S
  #
  KA <- rbind(P,R)
  KB <- rbind(Q,S)
  #
  # The F matrix is (KA,KB*N;0,N)
  F0 <- matrix(0,nrow=nrow(KA),ncol=nrow(R))
  #
  F1 <- cbind(KA,F0,KB%*%N)
  F2 <- cbind(matrix(rep(0,nrow(KA)*nrow(N)),nrow=ncol(N)),N)
  F  <- rbind(F1,F2)
  #
  # G is the 'impact' matrix for the shocks
  G <- rbind(KB,diag(ncol(N)))
  #
  return=list(F=F,G=G)
}

statespace.SDSGE <- function(obj){
  #
  if(obj$sol_type==1){ #gensys
    class(obj) <- "gensys"
    ssmats <- statespace(obj)
  }else if(obj$sol_type==2){ #uhlig
    class(obj) <- "uhlig"
    ssmats <- statespace(obj)
  }else{
    stop("Unrecognized solution method when forming state space model.\n",call.=FALSE)
  }
  #
  return=list(F=ssmats$F,G=ssmats$G)
}