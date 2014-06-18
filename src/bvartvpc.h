/*################################################################################
  ##
  ##   R package BMR by Keith O'Hara Copyright (C) 2011, 2012, 2013, 2014
  ##   This file is part of the R package BMR.
  ##
  ##   The R package BMR is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 3 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package BMR is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/
#ifndef _bvartvpc_H
#define _bvartvpc_H
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <RcppArmadillo.h>
RcppExport SEXP tsprior( SEXP , SEXP , SEXP , SEXP , SEXP );
RcppExport SEXP BVARTVPReps( SEXP , SEXP , SEXP , SEXP , SEXP , SEXP , 
                             SEXP , SEXP , SEXP , SEXP , SEXP ,
                             SEXP , SEXP , SEXP , SEXP , SEXP ,
                             SEXP , SEXP );
RcppExport SEXP BVARTVPIRFs( SEXP , SEXP , SEXP , SEXP , 
                             SEXP , SEXP );
#endif