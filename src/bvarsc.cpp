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

#include "bvarsc.h"
using namespace Rcpp;

SEXP SBVARReps( SEXP mX, SEXP mY, SEXP md, SEXP mdX, SEXP myd, SEXP mZd, SEXP mPsiPr, 
                SEXP minvPsiVPr, SEXP mBPr, SEXP mBeta, SEXP minvBVPr, SEXP mSigma, 
                SEXP mSigmaML, SEXP mgamma, SEXP mTp, SEXP mM, SEXP mp, SEXP mburnin, 
                SEXP mkreps )
{
  arma::mat Y = as<arma::mat>(mY);
  arma::mat X = as<arma::mat>(mX);
  arma::mat d = as<arma::mat>(md);
  arma::mat dX = as<arma::mat>(mdX);
  arma::mat yd = as<arma::mat>(myd);
  arma::mat Zd = as<arma::mat>(mZd);
  //
  arma::mat PsiPr = as<arma::mat>(mPsiPr);
  arma::mat invPsiVPr = as<arma::mat>(minvPsiVPr);
  arma::mat BPr = as<arma::mat>(mBPr);
  arma::mat BetaKT = as<arma::mat>(mBeta);
  arma::mat invBVPr = as<arma::mat>(minvBVPr);
  arma::mat Sigma = as<arma::mat>(mSigma);
  arma::mat SigmaML = as<arma::mat>(mSigmaML);
  int gamma = as<int>(mgamma);
  int Tp = as<int>(mTp);
  int M = as<int>(mM);
  int p = as<int>(mp);
  int burnin = as<int>(mburnin);
  int kreps = as<int>(mkreps);
  //
  int q = 1;
  //
  arma::cube PsiDraws(1, M, kreps);
  PsiDraws.zeros();
  arma::cube BetaDraws(p*M, M, kreps);
  BetaDraws.zeros();
  arma::cube SigmaDraws(M, M, kreps);
  SigmaDraws.zeros();
  //
  int i,j;
  //
  arma::mat KPsiRand = arma::randn(M*q);
  arma::mat KBetaRand = arma::randn(M*M*p);
  arma::mat KSigmaRand = arma::randn(Sigma.n_rows, Tp + p*M + gamma);
  //
  arma::mat D(Tp,p+1);
  D.zeros();
  D(arma::span(),arma::span(0,0)) = d;
  D(arma::span(),arma::span(1,p)) = dX;
  //
  arma::mat UFixed(M + p*M*q, M*q);
  UFixed.zeros();
  UFixed(arma::span(0,M-1),arma::span()) = arma::eye(M,M);
  arma::mat U = UFixed;
  //
  // Psi
  //
  arma::mat xi = Y - X*BetaKT;
  U = UFixed;
  arma::mat BetaT = trans(BetaKT);
  for(j=1;j<=p;j++){
    U(arma::span((M*q*j),M*q*(j+1)-1),arma::span()) = BetaT(arma::span(),arma::span((M*(j-1)),M*j-1));
  }
  //
  arma::mat invSigma = inv_sympd(Sigma);
  arma::mat invPsiVPt = (trans(U)*kron(trans(D)*D,invSigma)*U) + invPsiVPr;
  arma::mat StackTemp = invSigma*trans(xi)*D;
  arma::colvec StackedZD(StackTemp.begin(),StackTemp.size(),false);
  arma::mat PsiVPt = inv_sympd(invPsiVPt);
  arma::mat vecPsiPt = PsiVPt*(trans(U)*StackedZD + invPsiVPr*trans(PsiPr));
  KPsiRand = arma::randn(M*q);
  arma::mat vecPsi = vecPsiPt + trans(arma::chol(PsiVPt))*KPsiRand;
  arma::mat Psi(vecPsi.begin(),1,M,false);
  //
  // Beta
  yd = Y - d*Psi;
  Zd = X + dX*kron(arma::eye(p,p),Psi);
  arma::mat BetaVPt = inv_sympd(trans(Zd)*Zd + invBVPr);
  arma::mat BetaPt = BetaVPt*(trans(Zd)*yd + invBVPr*BPr);
  KBetaRand = arma::randn(M*M*p);
  arma::colvec vecBetaT(BetaPt.begin(),BetaPt.size(),false); 
  arma::mat vecBeta = vecBetaT + trans(arma::chol(kron(Sigma,BetaVPt)))*KBetaRand;
  arma::mat Beta(vecBeta.begin(),M*p,M);
  //
  // Sigma
  arma::mat Epsilon = yd - Zd*Beta;
  arma::mat B = trans(Beta - BPr)*invBVPr*(Beta - BPr);
  Sigma = trans(Epsilon)*Epsilon + SigmaML + B;
  KSigmaRand = arma::randn(Sigma.n_rows,Tp + p*M + gamma);
  Sigma = trans(arma::chol(inv_sympd(Sigma)))*KSigmaRand;
  Sigma = inv_sympd(Sigma*trans(Sigma));
  //
  //
  //
  for (i=1;i<=burnin;i++) {
    xi = Y - X*Beta;
    U = UFixed;
    arma::mat BetaT = trans(Beta);
    for(j=1;j<=p;j++){
      U(arma::span((M*q*j),M*q*(j+1)-1),arma::span()) = BetaT(arma::span(),arma::span((M*(j-1)),M*j-1));
    }
    //
    invSigma = inv_sympd(Sigma);
    invPsiVPt = (trans(U)*kron(trans(D)*D,invSigma)*U) + invPsiVPr;
    StackTemp = invSigma*trans(xi)*D;
    arma::colvec StackedZD(StackTemp.begin(),StackTemp.size(),false);
    PsiVPt = inv_sympd(invPsiVPt);
    vecPsiPt = PsiVPt*(trans(U)*StackedZD + invPsiVPr*trans(PsiPr));
    KPsiRand = arma::randn(M*q);
    vecPsi = vecPsiPt + trans(arma::chol(PsiVPt))*KPsiRand;
    arma::mat Psi(vecPsi.begin(),1,M,false);
    //
    // Beta
    yd = Y - d*Psi;
    Zd = X + dX*kron(arma::eye(p,p),Psi);
    BetaVPt = inv_sympd(trans(Zd)*Zd + invBVPr);
    BetaPt = BetaVPt*(trans(Zd)*yd + invBVPr*BPr);
    KBetaRand = arma::randn(M*M*p);
    arma::colvec vecBetaT(BetaPt.begin(),BetaPt.size(),false); 
    vecBeta = vecBetaT + trans(arma::chol(kron(Sigma,BetaVPt)))*KBetaRand;
    arma::mat Beta(vecBeta.begin(),M*p,M);
    //
    // Sigma
    Epsilon = yd - Zd*Beta;
    B = trans(Beta - BPr)*invBVPr*(Beta - BPr);
    Sigma = trans(Epsilon)*Epsilon + SigmaML + B;
    KSigmaRand = arma::randn(Sigma.n_rows,Tp + p*M + gamma);
    Sigma = trans(arma::chol(inv_sympd(Sigma)))*KSigmaRand;
    Sigma = inv_sympd(Sigma*trans(Sigma));
  }
  //
  for (i=1;i<=kreps;i++) {
    xi = Y - X*Beta;
    U = UFixed;
    arma::mat BetaT = trans(Beta);
    for(j=1;j<=p;j++){
      U(arma::span((M*q*j),M*q*(j+1)-1),arma::span()) = BetaT(arma::span(),arma::span((M*(j-1)),M*j-1));
    }
    //
    invSigma = inv_sympd(Sigma);
    invPsiVPt = (trans(U)*kron(trans(D)*D,invSigma)*U) + invPsiVPr;
    StackTemp = invSigma*trans(xi)*D;
    arma::colvec StackedZD(StackTemp.begin(),StackTemp.size(),false);
    PsiVPt = inv_sympd(invPsiVPt);
    vecPsiPt = PsiVPt*(trans(U)*StackedZD + invPsiVPr*trans(PsiPr));
    KPsiRand = arma::randn(M*q);
    vecPsi = vecPsiPt + trans(arma::chol(PsiVPt))*KPsiRand;
    arma::mat Psi(vecPsi.begin(),1,M,false);
    //
    // Beta
    yd = Y - d*Psi;
    Zd = X + dX*kron(arma::eye(p,p),Psi);
    BetaVPt = inv_sympd(trans(Zd)*Zd + invBVPr);
    BetaPt = BetaVPt*(trans(Zd)*yd + invBVPr*BPr);
    KBetaRand = arma::randn(M*M*p);
    arma::colvec vecBetaT(BetaPt.begin(),BetaPt.size(),false); 
    vecBeta = vecBetaT + trans(arma::chol(kron(Sigma,BetaVPt)))*KBetaRand;
    arma::mat Beta(vecBeta.begin(),M*p,M);
    //
    // Sigma
    Epsilon = yd - Zd*Beta;
    B = trans(Beta - BPr)*invBVPr*(Beta - BPr);
    Sigma = trans(Epsilon)*Epsilon + SigmaML + B;
    KSigmaRand = arma::randn(Sigma.n_rows,Tp + p*M + gamma);
    Sigma = trans(arma::chol(inv_sympd(Sigma)))*KSigmaRand;
    Sigma = inv_sympd(Sigma*trans(Sigma));
    //
    PsiDraws(arma::span(),arma::span(),arma::span(i-1,i-1)) = Psi;
    BetaDraws(arma::span(),arma::span(),arma::span(i-1,i-1)) = Beta;
    SigmaDraws(arma::span(),arma::span(),arma::span(i-1,i-1)) = Sigma;
  }
  //
  return Rcpp::List::create(Rcpp::Named("Psi") = PsiDraws,Rcpp::Named("Beta") = BetaDraws,Rcpp::Named("Sigma") = SigmaDraws);
}

SEXP SBVARIRFs( SEXP mM, SEXP mK, SEXP mkreps, SEXP mirfperiods, SEXP mBDs, SEXP mSDs )
{
  int M = as<int>(mM);
  int K = as<int>(mK);
  int kreps = as<int>(mkreps);
  int irfperiods = as<int>(mirfperiods);
  //
  NumericVector BDArray(mBDs);
  NumericVector SDArray(mSDs);
  arma::cube ImpStore(M, M, irfperiods*kreps);
  ImpStore.zeros();
  arma::cube BetaDraws(BDArray.begin(), K, M, kreps, false);
  arma::cube SigmaDraws(SDArray.begin(), M, M, kreps, false);
  //
  int i,j;
  //
  arma::mat BetaT(K,M);
  BetaT.zeros();
  arma::mat shock(M,M);
  shock.zeros();
  arma::mat ShockO(K-M,M);
  ShockO.zeros();
  arma::mat ShockB(K,M);
  ShockB.zeros();
  arma::mat Pow(M,M);
  Pow.zeros();
  //
  for (j=1;j<=kreps;j++) {
    BetaT = BetaDraws(arma::span(),arma::span(),arma::span(j-1,j-1));
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
      if (K > M){
        ShockB(arma::span(M,K-1),arma::span()) = ShockB(arma::span(0,K-M-1),arma::span());
      }
      ShockB(arma::span(0,M-1),arma::span()) = Pow;
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("ImpStore") = ImpStore);
}

SEXP bvarsforecast( SEXP mY0, SEXP mdX, SEXP mM, SEXP mp, SEXP mK, SEXP mruns, SEXP mperiods,
                    SEXP minclshocks, SEXP mPDs, SEXP mBDs, SEXP mSDs )
{
  int M = as<int>(mM);
  int p = as<int>(mp);
  int K = as<int>(mK);
  int runs = as<int>(mruns);
  int periods = as<int>(mperiods);
  int inclshocks = as<int>(minclshocks);
  //
  arma::mat Y0 = as<arma::mat>(mY0);
  arma::mat dX = as<arma::mat>(mdX);
  //
  NumericVector PDArray(mPDs);
  NumericVector BDArray(mBDs);
  NumericVector SDArray(mSDs);
  arma::cube PsiDraws(PDArray.begin(), 1, M, runs, false);
  arma::cube BetaDraws(BDArray.begin(), K, M, runs, false);
  arma::cube SigmaDraws(SDArray.begin(), M, M, runs, false);
  //
  arma::cube Forecasts(periods, M, runs);
  Forecasts.zeros();
  //
  arma::mat PT = PsiDraws(arma::span(),arma::span(),arma::span(0,0));
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
      PT = PsiDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
      BT = BetaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
      ST = SigmaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
      Shock = trans(arma::chol(ST));
      fY.zeros();
      kY = Y0;
      for (j=1;j<=periods;j++) {
        Krand = arma::randn(M);
        jshock = trans(Shock*Krand);
        //
        fY(arma::span(j-1,j-1),arma::span()) = kY*BT + PT - dX*(arma::kron(arma::eye(p,p),PT))*BT + jshock;
        if (K > M){
          kY(0,arma::span(M,K-1)) = kY(0,arma::span(0,K-M-1));
        }
        kY(0,arma::span(0,M-1)) = fY(arma::span(j-1,j-1),arma::span());
      }
      Forecasts(arma::span(),arma::span(),arma::span(i-1,i-1)) = fY;
    }
  }
  else{
    for (i=1;i<=runs;i++) {
      PT = PsiDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
      BT = BetaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
      ST = SigmaDraws(arma::span(),arma::span(),arma::span(i-1,i-1));
      fY.zeros();
      kY = Y0;
      for (j=1;j<=periods;j++) {
        fY(arma::span(j-1,j-1),arma::span()) = kY*BT + PT - dX*(arma::kron(arma::eye(p,p),PT))*BT;
        if (K > M){
          kY(0,arma::span(M,K-1)) = kY(0,arma::span(0,K-M-1));
        }
        kY(0,arma::span(0,M-1)) = fY(arma::span(j-1,j-1),arma::span());
      }
      Forecasts(arma::span(),arma::span(),arma::span(i-1,i-1)) = fY;
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("Forecasts") = Forecasts);
}