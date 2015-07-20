/*################################################################################
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
  ################################################################################*/
#ifndef _dsgec_H
#define _dsgec_H
#include <R.h>
//#include <Rmath.h>
#include <math.h>
#include <RcppArmadillo.h>
//#include "armadillo_qz"
RcppExport SEXP UhligCpp(SEXP mA, SEXP mB, SEXP mC, SEXP mD,
                         SEXP mF, SEXP mG, SEXP mH,
                         SEXP mJ, SEXP mK, SEXP mL, SEXP mM, SEXP mN,
                         SEXP mwhichEig);
RcppExport SEXP gensysCpp(SEXP mGamma0, SEXP mGamma1, SEXP mC, SEXP mPsi, SEXP mPi);
RcppExport SEXP DSGEKalman(SEXP dsgedata, SEXP mObserveMat, SEXP mObsCons,
                           SEXP mF, SEXP mG, SEXP mshocks,
                           SEXP mR, SEXP mmax_iter);
RcppExport SEXP DSGEKalmanFilt(SEXP dsgedata, SEXP mObserveMat,
                               SEXP mObsCons, SEXP mF, SEXP mG, SEXP mshocks,
                               SEXP mR, SEXP mmax_iter);
RcppExport SEXP DSGECR(SEXP dsgedata, SEXP mObserveMat,
                       SEXP mObsCons, SEXP mF, SEXP mG, SEXP mshocks,
                       SEXP mR, SEXP mmax_iter);
#endif