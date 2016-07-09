################################################################################
##
##   R package BMR by Keith O'Hara Copyright (C) 2011-2016
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

BVARM <- function(mydata,data_ext=NULL,coef_prior=NULL,constant=TRUE,p=4,n_draws=10000,VType=1,decay="H",HP1=0.5,HP2=0.5,HP3=1,HP4=2) UseMethod("BVARM")
BVARS <- function(mydata,psiprior=NULL,coefprior=NULL,p=4,irf.periods=20,keep=10000,burnin=1000,XiPsi=1,HP1=0.5,HP4=2,gamma=NULL) UseMethod("BVARS")
BVARW <- function(mydata,data_ext=NULL,coef_prior=NULL,constant=TRUE,p=4,n_draws=10000,n_burnin=1000,XiBeta=1,XiSigma=1,gamma=NULL) UseMethod("BVARW")
CVAR  <- function(mydata,data_ext,constant=TRUE,p=4,n_draws=10000) UseMethod("CVAR")
BVARTVP <- function(mydata,timelab=NULL,coefprior=NULL,tau=NULL,p=4,irf.periods=20,irf.points=NULL,keep=10000,burnin=5000,XiBeta=1,XiQ=0.01,gammaQ=NULL,XiSigma=1,gammaS=NULL) UseMethod("BVARTVP")

forecast <- function(obj, ...) UseMethod("forecast")
IRF <- function(obj, ...) UseMethod("IRF")
#
SDSGE <- function(mats, type=NULL) UseMethod("SDSGE")
gensys <- function(Gamma0,Gamma1,C,Psi,Pi) UseMethod("gensys")
uhlig <- function(A,B,C,D,F,G,H,J,K,L,M,N,whichEig=NULL) UseMethod("uhlig")

statespace <- function(obj) UseMethod("statespace")
DSGESim <- function(obj, ...) UseMethod("DSGESim")
# break up the following inputs into lists...
EDSGE <- function(dsgedata,chains=1,cores=1,ObserveMat,initialvals,partomats,priorform,priorpars,parbounds,parnames=NULL,optimMethod="Nelder-Mead",optimLower=NULL,optimUpper=NULL,optimControl=list(),DSGEIRFs=TRUE,irf.periods=20,scalepar=1,keep=50000,burnin=10000,tables=TRUE) UseMethod("EDSGE")
DSGEVAR <- function(dsgedata,chains=1,cores=1,lambda=Inf,p=2,constant=FALSE,ObserveMat,initialvals,partomats,priorform,priorpars,parbounds,parnames=NULL,optimMethod="Nelder-Mead",optimLower=NULL,optimUpper=NULL,optimControl=list(),IRFs=TRUE,irf.periods=20,scalepar=1,keep=50000,burnin=10000,tables=TRUE) UseMethod("DSGEVAR")

modecheck <- function(obj, ...) UseMethod("modecheck")
states <- function(obj, ...) UseMethod("states")
#
