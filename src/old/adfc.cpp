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

#include "adfc.h"
using namespace Rcpp;

SEXP ADFCheck( SEXP mY, SEXP mX, SEXP mp, SEXP mkT )
{
  arma::mat Y = as<arma::mat>(mY);
  arma::mat X = as<arma::mat>(mX);
  float kT = as<int>(mkT);
  int p = as<int>(mp);
  float K = 3;
  //
  //
  arma::mat Beta(1,1);
  Beta.zeros();
  arma::mat SigmaS(1,1);
  SigmaS.zeros();
  arma::mat SigmaS3 = SigmaS;
  float SigmaS2 = 1/(kT-1);
  arma::mat Epsilon = Y;
  Epsilon.zeros();
  //
  int i,j;
  arma::mat XT = X;
  //
  arma::mat LogLikelihood(1,1);
  LogLikelihood.zeros();
  arma::mat BIC(p,1);
  BIC.zeros();
  arma::uvec indices = arma::sort_index(BIC);
  arma::mat BICFinal(3,1);
  BICFinal.zeros();
  //
  arma::mat BICMin(1,1);
  BICMin.zeros();
  for(i=1;i<=3;i++){
    for(j=1;j<=p;j++){
      XT = X(arma::span(),arma::span(i-1,j+2));
      Beta = arma::solve(XT,Y);
      K = Beta.n_rows;
      Epsilon = Y - XT*Beta;
      SigmaS = (trans(Epsilon)*Epsilon);
      SigmaS2 = 1/(kT-K);
      SigmaS3 = SigmaS2*SigmaS;
      LogLikelihood = -((kT*(log(2*3.14159) + log(SigmaS3))) + inv(SigmaS3)*(SigmaS));
      LogLikelihood = 0.5*arma::as_scalar(LogLikelihood);
      BIC(arma::span(j-1,j-1),arma::span(0,0)) = -2*LogLikelihood + K*log(kT);
    }
    indices = arma::sort_index(BIC);
    BICMin = arma::as_scalar(indices(arma::span(0,0),arma::span(0,0)));
    BICFinal(arma::span(i-1,i-1),arma::span()) = BICMin;
  }
  //
  return Rcpp::List::create(Rcpp::Named("BIC") = BICFinal);
}