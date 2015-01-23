/*################################################################################
  ##
  ##   R package BMR by Keith O'Hara Copyright (C) 2011, 2012, 2013, 2014, 2015
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

#include "dsgevarc.h"
using namespace Rcpp;

SEXP DSGEVARPriorC( SEXP mdsgedata , SEXP mObserveMat , SEXP mObsCons , SEXP mF , SEXP mG , 
                    SEXP mN , SEXP mshocks , SEXP mR , SEXP mp , SEXP mMaxIter )
{
  try {
    arma::mat dsgedata = as<arma::mat>(mdsgedata);
    arma::mat ObserveMat = as<arma::mat>(mObserveMat);
    arma::mat ObsCons = as<arma::mat>(mObsCons);
    arma::mat F = as<arma::mat>(mF);
    arma::mat G = as<arma::mat>(mG);
    arma::mat N = as<arma::mat>(mN);
    arma::mat shocks = as<arma::mat>(mshocks);
    arma::mat R = as<arma::mat>(mR);
    //
    arma::mat OC2 = ObsCons*trans(ObsCons);
    //
    int p = as<int>(mp);
    //
    double MaxIter = as<int>(mMaxIter);
    //
    int i,j;
    //
    arma::mat Q(G.n_rows,G.n_rows);
    Q.zeros();
    Q(arma::span(G.n_rows-N.n_rows,G.n_rows-1),arma::span(G.n_rows-N.n_rows,G.n_rows-1)) = shocks;
    //
    arma::mat GQG = G*Q*trans(G);
    //
    arma::mat SigmaSS = arma::eye(F.n_cols,F.n_cols);
    //
    arma::mat SigmaSSOld = GQG;
    arma::mat SigmaSSNew = GQG;
    arma::mat IMat = F;
    //
    for(j=1;j<=MaxIter;j++){
      SigmaSSNew = SigmaSSOld + (IMat*SigmaSSOld*trans(IMat));
      IMat *= IMat;
      //
      SigmaSSOld = SigmaSSNew;
    }
    SigmaSS = SigmaSSNew;
    //
    arma::cube SigmaX(dsgedata.n_cols, dsgedata.n_cols, p+1);
    SigmaX.zeros();
    //
    SigmaX(arma::span(),arma::span(),arma::span(0,0)) = OC2 + trans(ObserveMat)*SigmaSS*ObserveMat + R;
    //
    arma::mat Fp = arma::eye(F.n_cols,F.n_cols);
    for(i=1;i<=p;i++){
      Fp *= F;
      SigmaX(arma::span(),arma::span(),arma::span(i,i)) = OC2 + trans(ObserveMat)*Fp*SigmaSS*ObserveMat;
    }
    //
    return Rcpp::List::create(Rcpp::Named("SigmaX") = SigmaX);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "BMR: DSGEVAR Prior C++ exception (unknown reason)" );
  }
  return R_NilValue;
}

SEXP DSGEVARLikelihood( SEXP mlogGPR, SEXP mXX, SEXP mGammaYY, SEXP mGammaXY , SEXP mGammaXX , 
                        SEXP mGammaBarYY , SEXP mGammaBarXY, SEXP mGammaBarXX , 
                        SEXP mlambda , SEXP mTp , SEXP mM , SEXP mp , SEXP mkcons )
{
  try {
    //
    float logGPR = as<double>(mlogGPR);
    float lambda = as<double>(mlambda);
    float Tp = as<int>(mTp);
    float M = as<int>(mM);
    float p = as<int>(mp);
    int kcons = as<int>(mkcons);
    //
    arma::mat XX = as<arma::mat>(mXX);
    //
    arma::mat GammaYY = as<arma::mat>(mGammaYY);
    arma::mat GammaXY = as<arma::mat>(mGammaXY);
    arma::mat GammaXX = as<arma::mat>(mGammaXX);
    //
    arma::mat GammaBarYY = as<arma::mat>(mGammaBarYY);
    arma::mat GammaBarXY = as<arma::mat>(mGammaBarXY);
    arma::mat GammaBarXX = as<arma::mat>(mGammaBarXX);
    //
    float lambdaT = Tp*lambda;
    float invlambda = 1/lambda;
    float M2 = M/2;
    //
    arma::mat invGammaXX = arma::inv_sympd(GammaXX);
    arma::mat invGammaBarXX = arma::inv_sympd(GammaBarXX);
    arma::mat SigmaEpsilon = GammaYY - GammaXY*invGammaXX*trans(GammaXY);
    arma::mat SigmaHatEpsilon = GammaBarYY - trans(GammaBarXY)*invGammaBarXX*GammaBarXY;
    #
    arma::mat invlambdaXX = invlambda*XX;
    double Term1 = -M2*log(arma::det(GammaXX + invlambdaXX)) + M2*log(arma::det(GammaXX));
    arma::mat invlambdaSigma = (1 + invlambda)*SigmaHatEpsilon;
    double Term2 = -((Tp + lambdaT - M*p - kcons)/2)*log(arma::det( invlambdaSigma ));
    double Term3 = ((lambdaT - M*p - kcons)/2)*log(arma::det(SigmaEpsilon));
    double logLikelihood = Term1 + Term2 + Term3 + logGPR - (((M*Tp)/2)*log(lambdaT*arma::datum::pi));
    //
    return Rcpp::List::create(Rcpp::Named("logLikelihood") = logLikelihood);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "BMR: DSGEVAR Likelihood C++ exception (unknown reason)" );
  }
  return R_NilValue;
}

SEXP DSGEVARLikelihoodInf( SEXP mYY , SEXP mXY , SEXP mXX , 
                           SEXP mGammaYY , SEXP mGammaXY , SEXP mGammaXX , 
                           SEXP mTp , SEXP mM , SEXP mp )
{
  try {
    //
    float Tp = as<int>(mTp);
    float M = as<int>(mM);
    float p = as<int>(mp);
    //
    arma::mat YY = as<arma::mat>(mYY);
    arma::mat XY = as<arma::mat>(mXY);
    arma::mat XX = as<arma::mat>(mXX);
    //
    arma::mat GammaYY = as<arma::mat>(mGammaYY);
    arma::mat GammaXY = as<arma::mat>(mGammaXY);
    arma::mat GammaXX = as<arma::mat>(mGammaXX);
    //
    arma::mat invGammaXX = arma::inv_sympd(GammaXX);
    arma::mat SigmaEpsilon = GammaYY - GammaXY*invGammaXX*trans(GammaXY);
    arma::mat Beta = invGammaXX*trans(GammaXY);
    arma::mat SigmaBarEpsilon = YY + trans(Beta)*XX*Beta - trans(Beta)*XY - trans(XY)*Beta;
    #
    double logLikelihood = - (((M*Tp)/2)*log(2*arma::datum::pi)) - ((Tp/2)*log(arma::det(SigmaEpsilon))) - ((Tp/2)*arma::trace(arma::inv_sympd(SigmaEpsilon)*SigmaBarEpsilon));
    //
    return Rcpp::List::create(Rcpp::Named("logLikelihood") = logLikelihood);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "BMR: DSGEVAR Likelihood C++ exception (unknown reason)" );
  }
  return R_NilValue;
}

SEXP DSGEVARReps( SEXP mGammaBarYY , SEXP mGammaBarXY , SEXP mGammaBarXX , SEXP mGXX , 
                  SEXP mGHXX , SEXP mlambda , SEXP mkreps , SEXP mTp , SEXP mM , SEXP mp ,
                  SEXP mkcons )
{
  try {
    //
    float lambda = as<double>(mlambda);
    int kreps = as<int>(mkreps);
    float Tp = as<int>(mTp);
    float M = as<int>(mM);
    float p = as<int>(mp);
    int kcons = as<int>(mkcons);
    //
    NumericVector GGammaBarYY(mGammaBarYY);
    NumericVector GGammaBarXY(mGammaBarXY);
    NumericVector GGammaBarXX(mGammaBarXX);
    NumericVector GGXX(mGXX);
    //
    arma::cube GammaBarYY(GGammaBarYY.begin(), M, M, kreps, false);
    arma::cube GammaBarXY(GGammaBarXY.begin(), M*p + kcons, M, kreps, false);
    arma::cube GammaBarXX(GGammaBarXX.begin(), M*p + kcons, M*p + kcons, kreps, false);
    arma::cube GXX(GGXX.begin(), M*p + kcons, M*p + kcons, kreps, false);
    //
    arma::mat GHXX = as<arma::mat>(mGHXX);
    //
    arma::cube BetaDraws(M*p + kcons, M, kreps);
    arma::cube SigmaDraws(M, M, kreps);
    //
    int j;
    //
    arma::mat GBarYY = GammaBarYY(arma::span(),arma::span(),arma::span(0,0));
    arma::mat GBarXY = GammaBarXY(arma::span(),arma::span(),arma::span(0,0));
    arma::mat GBarXX = GammaBarXX(arma::span(),arma::span(),arma::span(0,0));
    //
    arma::mat invGBarXX = arma::inv_sympd(GBarXX);
    arma::mat BetaHat = invGBarXX*GBarXY;
    arma::mat SigmaEpsilon = GBarYY - trans(GBarXY)*invGBarXX*GBarXY;
    //
    float S2 = (1+lambda)*Tp;
    arma::mat SigmaDraw = S2*SigmaEpsilon;
    int S3 = (1+lambda)*Tp - M*p - kcons;
    arma::mat KSigmaRand = arma::randn(SigmaDraw.n_rows,S3);
    arma::mat A = trans(arma::chol(arma::inv_sympd(SigmaDraw)))*KSigmaRand;
    SigmaDraw = arma::inv_sympd(A*trans(A));
    //
    float B1 = 1/lambda;
    arma::mat GammaXX = GXX(arma::span(),arma::span(),arma::span(0,0));
    arma::mat BetaCov = arma::inv_sympd(GammaXX + B1*GHXX);
    arma::vec KBetaRand = arma::randn(M*(M*p+kcons));
    arma::colvec VecBetaHat(BetaHat.begin(),BetaHat.size(),false);
    float B2 = 1/(Tp*lambda);
    arma::mat VecBetaDraw = VecBetaHat + trans(arma::chol(B2*arma::kron(SigmaEpsilon,BetaCov)))*KBetaRand;
    arma::mat BetaDraw(VecBetaDraw.begin(),M*p+kcons,M,false);
    //
    BetaDraws(arma::span(),arma::span(),arma::span(0,0)) = BetaDraw;
    SigmaDraws(arma::span(),arma::span(),arma::span(0,0)) = SigmaDraw;
    //
    for (j=2;j<=kreps;j++) {
      GBarYY = GammaBarYY(arma::span(),arma::span(),arma::span(j-1,j-1));
      GBarXY = GammaBarXY(arma::span(),arma::span(),arma::span(j-1,j-1));
      GBarXX = GammaBarXX(arma::span(),arma::span(),arma::span(j-1,j-1));
      //
      invGBarXX = arma::inv_sympd(GBarXX);
      BetaHat = invGBarXX*GBarXY;
      SigmaEpsilon = GBarYY - trans(GBarXY)*invGBarXX*GBarXY;
      //
      SigmaDraw = S2*SigmaEpsilon;
      KSigmaRand = arma::randn(SigmaDraw.n_rows,S3);
      A = trans(arma::chol(arma::inv(SigmaDraw)))*KSigmaRand;
      SigmaDraw = arma::inv_sympd(A*trans(A));
      //
      GammaXX = GXX(arma::span(),arma::span(),arma::span(j-1,j-1));
      BetaCov = arma::inv_sympd(GammaXX + B1*GHXX);
      KBetaRand = arma::randn(M*(M*p+kcons));
      arma::colvec VecBetaHat(BetaHat.begin(),BetaHat.size(),false);
      VecBetaDraw = VecBetaHat + trans(arma::chol(B2*arma::kron(SigmaEpsilon,BetaCov)))*KBetaRand;
      arma::mat BetaDraw(VecBetaDraw.begin(),M*p+kcons,M,false);
      //
      BetaDraws(arma::span(),arma::span(),arma::span(j-1,j-1)) = BetaDraw;
      SigmaDraws(arma::span(),arma::span(),arma::span(j-1,j-1)) = SigmaDraw;
    }
    //
    return Rcpp::List::create(Rcpp::Named("Beta") = BetaDraws,Rcpp::Named("Sigma") = SigmaDraws);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "BMR: DSGEVAR Gibbs C++ exception (unknown reason)" );
  }
  return R_NilValue;
}

SEXP DSGEVARRepsInf( SEXP mGammaBarYY , SEXP mGammaBarXY , SEXP mGammaBarXX , SEXP mlambda , 
                     SEXP mkreps , SEXP mTp , SEXP mM , SEXP mp , SEXP mkcons )
{
  try {
    //
    float lambda = as<double>(mlambda);
    int kreps = as<int>(mkreps);
    float Tp = as<int>(mTp);
    float M = as<int>(mM);
    float p = as<int>(mp);
    int kcons = as<int>(mkcons);
    //
    NumericVector GGammaBarYY(mGammaBarYY);
    NumericVector GGammaBarXY(mGammaBarXY);
    NumericVector GGammaBarXX(mGammaBarXX);
    //
    arma::cube GammaBarYY(GGammaBarYY.begin(), M, M, kreps, false);
    arma::cube GammaBarXY(GGammaBarXY.begin(), M*p + kcons, M, kreps, false);
    arma::cube GammaBarXX(GGammaBarXX.begin(), M*p + kcons, M*p + kcons, kreps, false);
    //
    arma::cube BetaDraws(M*p + kcons, M, kreps);
    arma::cube SigmaDraws(M, M, kreps);
    //
    int j;
    //
    arma::mat GBarYY = GammaBarYY(arma::span(),arma::span(),arma::span(0,0));
    arma::mat GBarXY = GammaBarXY(arma::span(),arma::span(),arma::span(0,0));
    arma::mat GBarXX = GammaBarXX(arma::span(),arma::span(),arma::span(0,0));
    //
    arma::mat invGBarXX = arma::inv_sympd(GBarXX);
    arma::mat BetaHat = invGBarXX*GBarXY;
    arma::mat SigmaEpsilon = GBarYY - trans(GBarXY)*invGBarXX*GBarXY;
    //
    BetaDraws(arma::span(),arma::span(),arma::span(0,0)) = BetaHat;
    SigmaDraws(arma::span(),arma::span(),arma::span(0,0)) = SigmaEpsilon;
    //
    for (j=2;j<=kreps;j++) {
      GBarYY = GammaBarYY(arma::span(),arma::span(),arma::span(j-1,j-1));
      GBarXY = GammaBarXY(arma::span(),arma::span(),arma::span(j-1,j-1));
      GBarXX = GammaBarXX(arma::span(),arma::span(),arma::span(j-1,j-1));
      //
      invGBarXX = arma::inv_sympd(GBarXX);
      BetaHat = invGBarXX*GBarXY;
      SigmaEpsilon = GBarYY - trans(GBarXY)*invGBarXX*GBarXY;
      //
      BetaDraws(arma::span(),arma::span(),arma::span(j-1,j-1)) = BetaHat;
      SigmaDraws(arma::span(),arma::span(),arma::span(j-1,j-1)) = SigmaEpsilon;
    }
    //
    return Rcpp::List::create(Rcpp::Named("Beta") = BetaDraws,Rcpp::Named("Sigma") = SigmaDraws);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "BMR: DSGEVAR Gibbs C++ exception (unknown reason)" );
  }
  return R_NilValue;
}

SEXP DSGEVARIRFs( SEXP mM , SEXP mK , SEXP mkcons , SEXP mkreps , SEXP mirfperiods , 
                  SEXP mBDs , SEXP mSDs , SEXP mIRF )
{
  try {
    int M = as<int>(mM);
    int K = as<int>(mK);
    int kcons = as<int>(mkcons);
    int kreps = as<int>(mkreps);
    int irfperiods = as<int>(mirfperiods);
    //
    NumericVector IRFArray(mIRF);
    NumericVector BDArray(mBDs);
    NumericVector SDArray(mSDs);
    arma::cube ImpStore(IRFArray.begin(), M, M, irfperiods*kreps, false);
    arma::cube BetaDraws(BDArray.begin(), K, M, kreps, false);
    arma::cube SigmaDraws(SDArray.begin(), M, M, kreps, false);
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
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "BMR: DSGEVAR IRF C++ exception (unknown reason)" );
  }
  return R_NilValue;
}

SEXP dsgevarforecast( SEXP mY0 , SEXP mM , SEXP mp , SEXP mK , SEXP mkcons ,
                      SEXP mruns , SEXP mperiods , SEXP minclshocks , SEXP mBDs , SEXP mSDs )
{
  try {
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
    NumericVector SDArray(mSDs);
    arma::cube BetaDraws(BDArray.begin(), K, M, runs, false);
    arma::cube SigmaDraws(SDArray.begin(), M, M, runs, false);
    //
    arma::cube Forecasts(periods, M, runs);
    Forecasts.zeros();
    //
    arma::mat BT = BetaDraws(arma::span(),arma::span(),arma::span(0,0));
    arma::mat ST = SigmaDraws(arma::span(),arma::span(),arma::span(0,0));
    arma::mat Shock = trans(arma::chol(ST));
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
        ST = SigmaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
        Shock = trans(arma::chol(ST));
        fY.zeros();
        kY = Y0;
        for (j=1;j<=periods;j++) {
          Krand = arma::randn(M);
          jshock = trans(Shock*Krand);
          //
          fY(arma::span(j-1,j-1),arma::span()) = kY*BT + jshock;
          if (K > M+kcons){
            kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1)); 
          }
          kY(0,arma::span(kcons,M-1+kcons)) = fY(arma::span(j-1,j-1),arma::span());
        }
        Forecasts(arma::span(),arma::span(),arma::span(i-1,i-1)) = fY;
      }
    }
    else{
      for (i=1;i<=runs;i++) {
        BT = BetaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
        ST = SigmaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
        fY.zeros();
        kY = Y0;
        for (j=1;j<=periods;j++) {
          fY(arma::span(j-1,j-1),arma::span()) = kY*BT;
          if (K > M+kcons){
            kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1)); 
          }
          kY(0,arma::span(kcons,M-1+kcons)) = fY(arma::span(j-1,j-1),arma::span());
        }
        Forecasts(arma::span(),arma::span(),arma::span(i-1,i-1)) = fY;
      }
    }
    //
    return Rcpp::List::create(Rcpp::Named("Forecasts") = Forecasts);
  } catch( std::exception &ex ) {
    forward_exception_to_r( ex );
  } catch(...) {
    ::Rf_error( "BMR: DSGEVAR Forecast C++ exception (unknown reason)" );
  }
  return R_NilValue;
}