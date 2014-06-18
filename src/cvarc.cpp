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

#include "cvarc.h"
using namespace Rcpp;

SEXP CVARReps( SEXP mBeta, SEXP mSigma, SEXP mX, SEXP mY, SEXP mTp, SEXP mM, 
              SEXP mK, SEXP mkcons, SEXP mboot, SEXP mRDs )
{
  arma::mat Y = as<arma::mat>(mY);
  arma::mat X = as<arma::mat>(mX);
  arma::mat Beta = as<arma::mat>(mBeta);
  arma::mat Sigma = as<arma::mat>(mSigma);
  float Tp = as<int>(mTp);
  int M = as<int>(mM);
  float K = as<int>(mK);
  int kcons = as<int>(mkcons);
  int boot = as<int>(mboot);
  //
  arma::cube BetaDraws(K, M, boot);
  BetaDraws.zeros();
  arma::cube SigmaDraws(M, M, boot);
  SigmaDraws.zeros();
  NumericVector RDArray(mRDs);
  arma::cube ResidDraws(RDArray.begin(), Tp, M, boot, false);
  //
  arma::mat BetaS = Beta;
  arma::mat SigmaS = Sigma;
  //float SigmaS2 = Tp/((Tp-K)*(Tp-K));
  float SigmaS2 = 1/(Tp-K);
  arma::mat Epsilon = Y;
  Epsilon.zeros();
  //
  int i,j;
  arma::mat YT(Tp,M);
  YT.zeros();
  arma::mat XT(Tp,K);
  XT.zeros();
  arma::mat XTM(1,K);
  XTM.zeros();
  XTM(0,0) = 1;
  arma::mat RT = YT;
  //
  for (i=1;i<=boot;i++) {
    RT = ResidDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
    XTM = X(arma::span(0,0),arma::span());
    for (j=1;j<=Tp;j++) {
      XT(arma::span(j-1,j-1),arma::span()) = XTM;
      YT(arma::span(j-1,j-1),arma::span()) = XTM*Beta + RT(arma::span(j-1,j-1),arma::span());
      XTM(0,arma::span(M+kcons,K-1)) = XTM(0,arma::span(kcons,K-M-1));
      XTM(0,arma::span(kcons,M-1+kcons)) = YT(arma::span(j-1,j-1),arma::span());
    }
    BetaS = arma::solve(XT,YT);
    Epsilon = YT - XT*BetaS;
    SigmaS = (trans(Epsilon)*Epsilon);
    SigmaS = SigmaS2*SigmaS;
    BetaDraws(arma::span(),arma::span(),arma::span(i-1,i-1)) = BetaS;
    SigmaDraws(arma::span(),arma::span(),arma::span(i-1,i-1)) = SigmaS;
  }
  //
  return Rcpp::List::create(Rcpp::Named("Beta") = BetaDraws,Rcpp::Named("Sigma") = SigmaDraws);
}

SEXP CVARIRFs( SEXP mM, SEXP mK, SEXP mkcons, SEXP mkreps, SEXP mirfperiods, 
              SEXP mBDs, SEXP mSDs )
{
  int M = as<int>(mM);
  int K = as<int>(mK);
  int kcons = as<int>(mkcons);
  int kreps = as<int>(mkreps);
  int irfperiods = as<int>(mirfperiods);
  //
  NumericVector BDArray(mBDs);
  NumericVector SDArray(mSDs);
  arma::cube BetaDraws(BDArray.begin(), K, M, kreps, false);
  arma::cube SigmaDraws(SDArray.begin(), M, M, kreps, false);
  arma::cube ImpStore(M, M, irfperiods*kreps);
  ImpStore.zeros();
  //
  int i,j;
  //
  arma::mat BetaT(K-kcons,M);
  BetaT.zeros();
  arma::mat shock(M,M);
  shock.zeros();
  arma::mat ShockO(K-M-kcons,M);
  ShockO.zeros();
  arma::mat ShockB(K-kcons,M);
  ShockB.zeros();
  arma::mat Pow(M,M);
  Pow.zeros();
  //
  for (j=1;j<=kreps;j++) {
    BetaT = BetaDraws(arma::span(kcons,K-1),arma::span(),arma::span(j-1,j-1));
    BetaT = trans(BetaT);
    shock = SigmaDraws(arma::span(),arma::span(),arma::span(j-1,j-1));
    shock = chol(shock);
    shock = trans(shock);
    ShockB.zeros();
    ShockB(arma::span(0,M-1),arma::span()) = shock;
    ImpStore(arma::span(),arma::span(),arma::span((j-1)*irfperiods,(j-1)*irfperiods)) = shock;
    for (i=2;i<=irfperiods;i++){
      Pow = BetaT*ShockB;
      ImpStore(arma::span(),arma::span(),arma::span((j-1)*irfperiods + (i-1),(j-1)*irfperiods + (i-1))) = Pow;
      ShockB(arma::span(M,K-1-kcons),arma::span()) = ShockB(arma::span(0,K-M-1-kcons),arma::span());
      ShockB(arma::span(0,M-1),arma::span()) = Pow;
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("ImpStore") = ImpStore);
}

SEXP cvarforecast( SEXP mY0, SEXP mM, SEXP mp, SEXP mK, SEXP mkcons, SEXP mperiods, 
                   SEXP mci, SEXP mBeta, SEXP mSigma )
{
  int M = as<int>(mM);
  int p = as<int>(mp);
  int K = as<int>(mK);
  int kcons = as<int>(mkcons);
  int periods = as<int>(mperiods);
  float ci = as<float>(mci);
  //
  arma::mat Y0 = as<arma::mat>(mY0);
  //
  arma::mat Beta = as<arma::mat>(mBeta);
  arma::mat BetaNC = Beta(arma::span(0+kcons,K-1),arma::span());
  arma::mat Sigma = as<arma::mat>(mSigma);
  //
  arma::cube Forecasts(periods, M, 3);
  Forecasts.zeros();
  //
  arma::mat fY(periods,M);
  fY.zeros();
  arma::mat fYL(periods,M);
  fYL.zeros();
  arma::mat fYU(periods,M);
  fYU.zeros();
  arma::mat kY = Y0;
  //
  arma::vec SEY(M,1);
  SEY.zeros();
  //
  int i,j;
  //
  arma::mat Shocks(p*M,p*M);
  Shocks.zeros();
  arma::mat SEForecast(M,M);
  SEForecast.zeros();
  for (j=1;j<=periods;j++) {
    SEForecast = Sigma + trans(BetaNC)*Shocks*BetaNC;
    Shocks(arma::span(M,(M*p)-1),arma::span(M,(M*p)-1)) = Shocks(arma::span(0,(M*(p-1))-1),arma::span(0,(M*(p-1))-1));
    Shocks(arma::span(0,M-1),arma::span(0,M-1)) = Sigma;
    //
    SEY = SEForecast.diag();
    //
    fY(arma::span(j-1,j-1),arma::span()) = kY*Beta;
    fYL(arma::span(j-1,j-1),arma::span()) = kY*Beta - ci*trans(sqrt(SEY));
    fYU(arma::span(j-1,j-1),arma::span()) = kY*Beta + ci*trans(sqrt(SEY));
    //
    kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1)); 
    kY(0,arma::span(kcons,M-1+kcons)) = fY(arma::span(j-1,j-1),arma::span());
    //
  }
  // Point Forecast, then lower CI, then upper CI
  Forecasts(arma::span(),arma::span(),arma::span(0,0)) = fY;
  Forecasts(arma::span(),arma::span(),arma::span(1,1)) = fYL;
  Forecasts(arma::span(),arma::span(),arma::span(2,2)) = fYU;
  //
  return Rcpp::List::create(Rcpp::Named("Forecasts") = Forecasts,Rcpp::Named("Shocks") = Shocks);
}