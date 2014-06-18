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

#include "bvarmc.h"
using namespace Rcpp;

SEXP MBVARReps( SEXP mSigma, SEXP mZ, SEXP mY, SEXP maPr, SEXP mBVPr, SEXP mM, SEXP mK, 
                SEXP mburnin, SEXP mkreps )
{
  arma::mat Sigma = as<arma::mat>(mSigma);
  arma::mat Z = as<arma::mat>(mZ);
  arma::mat Y = as<arma::mat>(mY);
  arma::mat aPr = as<arma::mat>(maPr);
  arma::mat BVPr = as<arma::mat>(mBVPr);
  int M = as<int>(mM);
  int K = as<int>(mK);
  int burnin = as<int>(mburnin);
  int kreps = as<int>(mkreps);
  //
  arma::cube BetaDraws(K, M, kreps);
  BetaDraws.zeros();
  //
  int i,j;
  //
  arma::mat invBVPr = inv(BVPr);
  arma::mat BVPt = inv_sympd(invBVPr + trans(Z)*Sigma*Z);
  arma::colvec YStacked(Y.begin(), Y.size(), false);
  arma::mat aPt = BVPt*(invBVPr*aPr + trans(Z)*Sigma*YStacked);
  arma::vec Krand = arma::randn(M*K);
  arma::mat tCholBVPt = trans(arma::chol(BVPt));
  arma::vec alpha = aPt + tCholBVPt*Krand;
  arma::mat Beta(alpha.begin(),K,M,false);
  //
  for (i=1;i<=burnin;i++) {
    Krand = arma::randn(M*K);
    alpha = aPt + tCholBVPt*Krand;
  }
  //
  for (j=1;j<=kreps;j++) {
    Krand = arma::randn(M*K);
    alpha = aPt + tCholBVPt*Krand;
    arma::mat Beta(alpha.begin(),K,M,false);
    //
    BetaDraws(arma::span(),arma::span(),arma::span(j-1,j-1)) = Beta;
  }
  //
  return Rcpp::List::create(Rcpp::Named("Beta") = BetaDraws);
}

SEXP MBVARIRFs( SEXP mshock, SEXP mM, SEXP mK, SEXP mkcons, SEXP mkreps, 
                SEXP mirfperiods, SEXP mBDs )
{
  arma::mat shock = as<arma::mat>(mshock);
  int M = as<int>(mM);
  int K = as<int>(mK);
  int kcons = as<int>(mkcons);
  int kreps = as<int>(mkreps);
  int irfperiods = as<int>(mirfperiods);
  //
  NumericVector BDArray(mBDs);
  arma::cube ImpStore(M, M, irfperiods*kreps);
  ImpStore.zeros();
  arma::cube BetaDraws(BDArray.begin(), K, M, kreps, false);
  //
  int i,j;
  //
  shock = trans(shock);
  //
  arma::mat BetaT(K-kcons,M);
  BetaT.zeros();
  arma::mat ShockO(K-M-kcons,M);
  ShockO.zeros();
  arma::mat ShockB(K-kcons,M);
  ShockB.zeros();
  arma::mat Pow(M,M);
  Pow.zeros();
  //
  ShockB(arma::span(0,M-1),arma::span()) = shock;
  //
  for (j=1;j<=kreps;j++) {
    BetaT = BetaDraws(arma::span(kcons,K-1),arma::span(),arma::span(j-1,j-1));
    BetaT = trans(BetaT);
    ShockB.zeros();
    ShockB(arma::span(0,M-1),arma::span()) = shock;
    ImpStore(arma::span(),arma::span(),arma::span((j-1)*irfperiods,(j-1)*irfperiods)) = shock;
    for (i=2;i<=irfperiods;i++){
      Pow = BetaT*ShockB;
      ImpStore(arma::span(),arma::span(),arma::span((j-1)*irfperiods + (i-1),(j-1)*irfperiods + (i-1))) = Pow;
      if (K > M+kcons){
        ShockB(arma::span(M,K-1-kcons),arma::span()) = ShockB(arma::span(0,K-M-1-kcons),arma::span());
      }
      ShockB(arma::span(0,M-1),arma::span()) = Pow;
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("ImpStore") = ImpStore);
}

SEXP bvarmforecast( SEXP mY0, SEXP mM, SEXP mp, SEXP mK, SEXP mkcons, SEXP mruns, SEXP mperiods, 
                SEXP minclshocks, SEXP mBDs, SEXP mShock )
{
  int M = as<int>(mM);
  int p = as<int>(mp);
  int K = as<int>(mK);
  int kcons = as<int>(mkcons);
  int runs = as<int>(mruns);
  int periods = as<int>(mperiods);
  int inclshocks = as<int>(minclshocks);
  //
  arma::mat Y0 = as<arma::mat>(mY0);
  //
  NumericVector BDArray(mBDs);
  arma::cube BetaDraws(BDArray.begin(), K, M, runs, false);
  //
  arma::mat Shock = as<arma::mat>(mShock);
  //
  arma::cube Forecasts(periods, M, runs);
  Forecasts.zeros();
  //
  arma::mat BT = BetaDraws(arma::span(),arma::span(),arma::span(0,0));
  //
  arma::mat fY(periods,M);
  fY.zeros();
  arma::mat kY = Y0;
  arma::vec Krand = arma::randn(M);
  arma::mat jshock(1,M);
  jshock.zeros();
  //
  int i,j;
  //
  if(inclshocks > 0){
    for (i=1;i<=runs;i++) {
      BT = BetaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
      fY.zeros();
      kY = Y0;
      for (j=1;j<=periods;j++) {
        Krand = arma::randn(M);
        jshock = trans(Shock*Krand);
        //
        fY(arma::span(j-1,j-1),arma::span()) = kY*BT + jshock;
        kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1)); 
        kY(0,arma::span(kcons,M-1+kcons)) = fY(arma::span(j-1,j-1),arma::span());
      }
      Forecasts(arma::span(),arma::span(),arma::span(i-1,i-1)) = fY;
    }
  }
  else{
    for (i=1;i<=runs;i++) {
      BT = BetaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
      fY.zeros();
      kY = Y0;
      for (j=1;j<=periods;j++) {
        fY(arma::span(j-1,j-1),arma::span()) = kY*BT;
        kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1)); 
        kY(0,arma::span(kcons,M-1+kcons)) = fY(arma::span(j-1,j-1),arma::span());
      }
      Forecasts(arma::span(),arma::span(),arma::span(i-1,i-1)) = fY;
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("Forecasts") = Forecasts);
}