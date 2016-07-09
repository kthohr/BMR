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

CVAR.default <- function(mydata,data_ext,constant=TRUE,p=4,n_draws=10000)
{
    #
    err_check <- .cvar_errors(mydata,p)
    #
    cvar_obj <- .cvar_run(mydata,data_ext,constant,p,n_draws)
    #
    cvar_ret <- list(data=mydata,data_ext=data_ext,cvar_obj=cvar_obj)
    class(cvar_ret) <- "CVAR"
    #
    return(cvar_ret)
}

.cvar_errors <- function(mydata,p)
{
    #
    # basic sanity checks
    if(ncol(mydata)<2){
        stop("need more than 1 variable in your data.\n",call.=FALSE)
    }
    if(p<1){
        stop("need at least 1 lag; set p higher.\n",call.=FALSE)
    }
    if(p>nrow(mydata)){
        stop("need more data points than lags.\n",call.=FALSE)
    }
    for(i in 1:ncol(mydata)){
        if(sum(is.na(mydata[,i]))){
            stop("no missing observations allowed.\n",call.=FALSE)
        }
    }
    #
}

.cvar_run <- function(mydata,data_ext,constant,p,n_draws)
{
    cvar_obj <- new(R_cvar)
    data_raw <- as.matrix(mydata)
    #
    cvar_obj$cons_term = constant
    cvar_obj$p = p

    if (!is.null(data_ext)) {
        data_ext <- as.matrix(data_ext)
        cvar_obj$data(data_raw,data_ext)
    } else {
        cvar_obj$data(data_raw)
    }
    #
    cvar_obj$setup()
    #
    cvar_obj$boot(n_draws)
    #
    return(cvar_obj)
}