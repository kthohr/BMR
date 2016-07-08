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

# 07/07/2016

BVARW.default <- function(mydata,data_ext=NULL,coef_prior=NULL,constant=TRUE,p=4,n_draws=10000,n_burnin=1000,XiBeta=1,XiSigma=1,gamma=NULL)
{
    #
    err_check <- .bvarwerrors(mydata,coef_prior,constant,p,XiBeta,XiSigma,gamma)
    #
    bvar_obj <- .bvarw_run(mydata,data_ext,err_check$coef_prior,constant,p,n_draws,n_burnin,err_check$XiBeta,err_check$XiSigma,err_check$gamma)
    #
    bvarw_ret <- list(data=mydata,data_ext=data_ext,bvar_obj=bvar_obj)
    class(bvarw_ret) <- "BVARW"
    #
    return(bvarw_ret)
}

.bvarwerrors <- function(mydata,coef_prior,constant,p,XiBeta,XiSigma,gamma)
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

    for (i in 1:ncol(mydata)) {
        if (sum(is.na(mydata[,i]))) {
            stop("no missing observations allowed.\n",call.=FALSE)
        }
    }
    #
    # Errors around user-given priors. First, Beta.
    #
    #No priors specified for the coefficients, set equal to random walk (in levels)
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
    # Prior on the variance of alpha
    #
    if (class(XiBeta) == "numeric") {
        if (XiBeta <= 0) {
            stop("XiBeta must be greater than zero.\n",call.=FALSE)
        }
    } else if (class(XiBeta) == "matrix") {
        if (constant) {
            if (nrow(XiBeta) != ((ncol(mydata)*p+1)*ncol(mydata)) || ncol(XiBeta) != ((ncol(mydata)*p+1)*ncol(mydata))) {
                stop("you have selected a full matrix for XiBeta.\n","Therefore, XiBeta must be of dimensions ",((ncol(mydata)*p+1)*ncol(mydata))," x ",((ncol(mydata)*p+1)*ncol(mydata)),".\n",call.=FALSE)
            }
        } else {
            if (nrow(XiBeta) != ((ncol(mydata)*p)*ncol(mydata)) || ncol(XiBeta) != ((ncol(mydata)*p)*ncol(mydata))) {
                stop("you have selected a full matrix for XiBeta.\n","Therefore, XiBeta must be of dimensions ",((ncol(mydata)*p)*ncol(mydata))," x ",((ncol(mydata)*p)*ncol(mydata)),".\n",call.=FALSE)
            }
        }
    } else {
        if (constant) {
            stop("unrecognised form for XiBeta.\n","XiBeta must be a ",((ncol(mydata)*p+1)*ncol(mydata))," x ",((ncol(mydata)*p+1)*ncol(mydata))," matrix, or single numeric value.\n",call.=FALSE)
        } else {
            stop("unrecognised form for XiBeta.\n","XiBeta must be a ",((ncol(mydata)*p)*ncol(mydata))," x ",((ncol(mydata)*p)*ncol(mydata))," matrix, or single numeric value.\n",call.=FALSE)
        }
    }
    #
    # Now Sigma
    # 
    if (class(XiSigma) == "numeric") {
        if (XiSigma <= 0) {
            stop("XiSigma must be greater than zero.\n",call.=FALSE)
        }

        XiSigma <- diag(ncol(mydata)) * XiSigma
    } else if (class(XiSigma) == "matrix") {
        if (nrow(XiSigma) != (ncol(mydata)) || ncol(XiSigma) != (ncol(mydata))) {
            stop("you have selected a full matrix for XiSigma.\n","Therefore, XiSigma must be of dimensions ",(ncol(mydata))," x ",(ncol(mydata)),".\n",call.=FALSE)
        }
    } else {
        stop("unrecognised form for XiSigma.\n","XiSigma must be a ",(ncol(mydata))," x ",(ncol(mydata))," matrix, or single numeric value.\n",call.=FALSE)
    }
    #
    # Prior COV Degrees of Freedom
    if (class(gamma)=="NULL") {
        gamma <- ncol(mydata) + 1
    }
    
    if (gamma <= ncol(mydata)) {
        stop("with this data, the minimum value for gamma is ",(ncol(mydata)+1),".\n",call.=FALSE)
    }
    #
    #
    #
    return=list(coef_prior=coef_prior,XiBeta=XiBeta,XiSigma=XiSigma,gamma=gamma)
}

.bvarw_run <- function(mydata,data_ext,coef_prior,constant,p,n_draws,n_burnin,XiBeta,XiSigma,gamma)
{
    bvar_obj <- new(R_bvarw)
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
    bvar_obj$prior(coef_prior,XiBeta,XiSigma,gamma)
    #
    bvar_obj$gibbs(n_draws,n_burnin)
    #
    return(bvar_obj)
}

