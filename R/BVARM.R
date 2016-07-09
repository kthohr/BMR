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

BVARM.default <- function(mydata,data_ext=NULL,coef_prior=NULL,constant=TRUE,p=4,n_draws=10000,VType=1,decay="H",HP1=0.5,HP2=0.5,HP3=1,HP4=2)
{
    #
    err_check <- .bvarm_errors(mydata,p,coef_prior,constant,VType,decay,HP4)
    #
    bvar_obj <- .bvarm_run(mydata,data_ext,err_check$coef_prior,constant,p,n_draws,VType,decay,HP1,HP2,HP3,HP4)
    #
    bvarm_ret <- list(data=mydata,data_ext=data_ext,bvar_obj=bvar_obj)
    class(bvarm_ret) <- "BVARM"
    #
    return(bvarm_ret)
}

.bvarm_errors <- function(mydata,p,coef_prior,constant,VType,decay,HP4)
{
    #
    # basic sanity checks
    if (ncol(mydata) < 2) {
        stop("need more than 1 variable.\n",call.=FALSE)
    }

    if (p < 1) {
        stop("need at least 1 lag.\n",call.=FALSE)
    }

    if (p > nrow(mydata)) {
        stop("need more data points than lags.\n",call.=FALSE)
    }

    for(i in 1:ncol(mydata)){
        if(sum(is.na(mydata[,i]))){
            stop("no missing observations allowed.\n",call.=FALSE)
        }
    }
    #
    # Errors around prior mean of the coefficients
    #
    # if no prior is provided for the coefficients, set to a random walk (in levels)
    if (class(coef_prior)=="NULL") {
        coef_prior <- c(rep(1,ncol(mydata)))
    }

    if (class(coef_prior)=="matrix") {
        if (constant) {
            if (length(c(coef_prior)) != ((ncol(mydata)*p + 1)*ncol(mydata))) {
                stop("you have opted to give a full prior matrix on Beta. However,\n", "it is not of the appropriate dimensions, which are ",(ncol(mydata)*p + 1)," x ",ncol(mydata),".\n",call.=FALSE)
            }
        } else {
            if (length(c(coef_prior)) != ((ncol(mydata)*p)*ncol(mydata))) {
                stop("you have opted to give a full prior matrix on Beta. However,\n", "it is not of the appropriate dimensions, which are ",(ncol(mydata)*p)," x ",ncol(mydata),".\n",call.=FALSE)
            }
        }
    }

    if (class(coef_prior)=="numeric") {
        if(length(coef_prior) != ncol(mydata)){
            stop("you have opted to give a prior on the first-own lags in Beta. However,\n", "it is not of the appropriate dimensions; please provide a numeric vector of length ",ncol(mydata),".\n",call.=FALSE)
        }
    }
    #
    # Some possible errors around hyperparameters and such
    #
    if (VType != 1 && VType != 2) {
        stop("VType can only be 1 or 2.\n",call.=FALSE)
    }
    if (decay != "H" && decay != "G") {
        stop("decay must be H or G (and quotation marks both sides!).\n",call.=FALSE)
    }
    if (HP4 <= 0) {
        stop("HP4 must be greater than zero.\n",call.=FALSE)
    }
    #
    return=list(coef_prior=coef_prior)
}

.bvarm_run <- function(mydata,data_ext,coef_prior,constant,p,n_draws,VType,decay,HP1,HP2,HP3,HP4)
{
    bvar_obj <- new(R_bvarm)
    data_raw <- as.matrix(mydata)
    #
    bvar_obj$cons_term = constant
    bvar_obj$p = p

    if (!is.null(data_ext)) {
        data_ext <- as.matrix(data_ext)
        bvar_obj$data(data_raw,data_ext)
    } else {
        bvar_obj$data(data_raw)
    }
    #
    bvar_obj$var_type = VType

    if (decay=="G") {
        bvar_obj$decay_type = 1
    } else {
        bvar_obj$decay_type = 2
    }

    bvar_obj$hyper_pars = c(HP1,HP2,HP3,HP4)

    bvar_obj$prior(coef_prior)
    #
    bvar_obj$gibbs(n_draws)
    #
    return(bvar_obj)
}
