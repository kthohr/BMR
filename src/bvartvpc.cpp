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

#include "bvartvpc.h"
using namespace Rcpp;

SEXP tsprior( SEXP mY, SEXP mtau, SEXP mM, SEXP mK, SEXP mp )
{
  arma::mat Y = as<arma::mat>(mY);
  int tau = as<int>(mtau);
  int M = as<int>(mM);
  int K = as<int>(mK);
  int p = as<int>(mp);
  //
  arma::mat Z(M*tau,K);
  Z.zeros();
  //
  arma::mat BetaDenominator(K,K);
  BetaDenominator.zeros();
  arma::mat BetaNumerator(K,1);
  BetaNumerator.zeros();
  //
  int i,j,k;
  //
  //
  //
  arma::mat Ytau = trans(Y(arma::span(p,tau+p-1),arma::span()));
  //
  for(i = (p+1); i<= (tau+p);i++){
    k = i - p - 1;
    Z(arma::span(k*M,(k+1)*M-1),arma::span(0,M-1)) = arma::eye(M,M);
    for(j=1;j<=p;j++){
      Z(arma::span(k*M,(k+1)*M-1),arma::span(M+((M*M)*(j-1)),M-1+((M*M)*(j)))) = kron(arma::eye(M,M),Y(arma::span(i-1-j,i-1-j),arma::span()));
    }
  }
  //
  arma::mat zt(M,K);
  zt.zeros();
  arma::mat Beta(K,1);
  Beta.zeros();
  //
  for(i=1;i<=tau;i++){
    zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
    BetaDenominator = BetaDenominator + trans(zt)*zt;
    BetaNumerator = BetaNumerator + trans(zt)*Ytau(arma::span(),arma::span(i-1,i-1));
  }
  BetaDenominator = inv(BetaDenominator);
  Beta = BetaDenominator*BetaNumerator;
  //
  arma::mat SSE(M,M);
  SSE.zeros();
  //
  for(i=1;i<=tau;i++){
    zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
    SSE += (Ytau(arma::span(),arma::span(i-1,i-1)) - zt*Beta)*trans(Ytau(arma::span(),arma::span(i-1,i-1)) - zt*Beta);
  }
  //
  float tau2 = tau - (M*p + 1);
  arma::mat SPrOLS = SSE*(1/tau2);
  arma::mat SPrOLSinv = inv(SPrOLS);
  //
  arma::mat BVPrOLS(K,K);
  BVPrOLS.zeros();
  for(i=1;i<=tau;i++){
    zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
    BVPrOLS = BVPrOLS + trans(zt)*SPrOLSinv*zt;
  }
  BVPrOLS = inv(BVPrOLS);
  //
  return Rcpp::List::create(Rcpp::Named("BetaOLS") = Beta,Rcpp::Named("BVPrOLS") = BVPrOLS, Rcpp::Named("SigmaOLS") = SPrOLS);
}

SEXP BVARTVPReps( SEXP my, SEXP mZ, SEXP mM, SEXP mK, SEXP mkT, SEXP mkeep, 
                  SEXP mburnin, SEXP mB0Pr, SEXP mB0VPr, SEXP minvB0VPr, SEXP mQPr,
                  SEXP mQVPr, SEXP mSPr, SEXP mSVPr, SEXP mQDraw, SEXP mQChol,
                  SEXP mSDraw, SEXP minvSDraw )
{
  //
  arma::mat y = as<arma::mat>(my);
  arma::mat Z = as<arma::mat>(mZ);
  arma::mat B0Pr = as<arma::mat>(mB0Pr);
  arma::mat B0VPr = as<arma::mat>(mB0VPr);
  arma::mat invB0VPr = as<arma::mat>(minvB0VPr);
  arma::mat QPr = as<arma::mat>(mQPr);
  arma::mat SPr = as<arma::mat>(mSPr);
  arma::mat QDraw = as<arma::mat>(mQDraw);
  arma::mat QChol = as<arma::mat>(mQChol);
  arma::mat SigmaDraw = as<arma::mat>(mSDraw);
  arma::mat invSDraw = as<arma::mat>(minvSDraw);
  //
  int M = as<int>(mM);
  int K = as<int>(mK);
  int kT = as<int>(mkT);
  int keep = as<int>(mkeep);
  int burnin = as<int>(mburnin);
  int QVPr = as<int>(mQVPr);
  int SVPr = as<int>(mSVPr);
  //
  int i,j;
  int kreps;
  //
  arma::cube BetaDraws(K,(kT+1),keep);
  BetaDraws.zeros();
  arma::cube QDraws(K,K,keep);
  QDraws.zeros();
  arma::cube SDraws(M,M,keep);
  SDraws.zeros();
  //
  arma::mat BtDraw(K,(kT+1));
  BtDraw.zeros();
  arma::mat BetaTilde(K,(kT+1));
  BetaTilde.zeros();
  //
  arma::mat SChol = arma::chol(SigmaDraw);
  //
  arma::mat BVPt(K,K);
  BVPt.zeros();
  arma::mat BetaNumerator(K,1);
  BetaNumerator.zeros();
  arma::mat Zt(M,K);
  Zt.zeros();
  //
  for(i=1; i<=kT; i++){
    Zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
    BVPt += trans(Zt)*invSDraw*Zt;
    BetaNumerator += trans(Zt)*invSDraw*(y(arma::span(),arma::span(i-1,i-1)) - Zt*BetaTilde(arma::span(),arma::span(i-1,i-1)));
  }
  BVPt = inv_sympd(BVPt + invB0VPr);
  arma::mat KBRand = arma::randn(K,1);
  arma::mat B0Draw = BVPt*(invB0VPr*B0Pr + BetaNumerator) + trans(arma::chol(BVPt))*KBRand;
  //
  arma::mat Epsilon(M,kT);
  Epsilon.zeros();
  for(i=1; i<=kT; i++){
    Epsilon(arma::span(),arma::span(i-1,i-1)) = y(arma::span(),arma::span(i-1,i-1)) - Z(arma::span((i-1)*M,(i*M)-1),arma::span())*B0Draw;
  }
  //
  arma::mat ResidPlus((M+K)*kT,1);
  ResidPlus.zeros();
  //
  arma::mat KWRand1 = arma::randn(M,1);
  arma::mat KWRand2 = arma::randn(K,1);
  //
  for(i=1; i<=kT; i++){
    ResidPlus(arma::span((i-1)*(M+K),(i-1)*(M+K)+M-1),arma::span(0,0)) = trans(SChol)*KWRand1;
    ResidPlus(arma::span((i-1)*(M+K)+M,(i-1)*(M+K)+(M+K)-1),arma::span(0,0)) = trans(QChol)*KWRand2;
    //
    KWRand1 = arma::randn(M,1);
    KWRand2 = arma::randn(K,1);
  }
  //
  arma::mat BetaPlus(K,kT+1);
  BetaPlus.zeros();
  arma::mat YBetaPlus(M,kT);
  YBetaPlus.zeros();
  //
  for(i=1; i<=kT; i++) {
    Zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
    YBetaPlus(arma::span(),arma::span(i-1,i-1)) = Zt*BetaPlus(arma::span(),arma::span(i-1,i-1)) + ResidPlus(arma::span((i-1)*(M+K),(i-1)*(M+K)+M-1),arma::span());
    BetaPlus(arma::span(),arma::span(i,i)) = BetaPlus(arma::span(),arma::span(i-1,i-1)) + ResidPlus(arma::span((i-1)*(M+K) + M,(i*(M+K))-1),arma::span());
  }
  //
  arma::mat EpsilonYBetaPlus = Epsilon - YBetaPlus;
  //
  arma::mat LStore(K*kT,K);
  LStore.zeros();
  arma::mat FStore(M*kT,M);
  FStore.zeros();
  arma::mat BetaFiltered(K,kT+1);
  BetaFiltered.zeros();
  arma::mat KalmanResid(M,kT);
  KalmanResid.zeros();
  arma::mat StateCovFiltered(K,K);
  StateCovFiltered.zeros();
  //
  arma::mat rt(K,kT+1);
  rt.zeros();
  arma::mat BetaSmoothed(K,kT+1);
  BetaSmoothed.zeros();
  //
  Zt.zeros();
  //
  Zt = Z(arma::span(0,M-1),arma::span());
  KalmanResid(arma::span(),arma::span(0,0)) = EpsilonYBetaPlus(arma::span(),arma::span(0,0)) - Zt*BetaFiltered(arma::span(),arma::span(0,0));
  arma::mat StateCovPredicted = Zt*StateCovFiltered*trans(Zt) + SigmaDraw;
  arma::mat invStateCovPredicted = inv_sympd(StateCovPredicted);
  FStore(arma::span(0,M-1),arma::span()) = invStateCovPredicted;
  arma::mat KalmanGain = StateCovFiltered*trans(Zt)*invStateCovPredicted;
  arma::mat Lt = arma::eye(K,K) - KalmanGain*Zt;
  LStore(arma::span(0,K-1),arma::span()) = Lt;
  BetaFiltered(arma::span(),arma::span(1,1)) = BetaFiltered(arma::span(),arma::span(0,0)) + KalmanGain*KalmanResid(arma::span(),arma::span(0,0));
  StateCovFiltered = StateCovFiltered*trans(Lt) + QDraw;
  //
  for(i=2; i<=kT; i++) {
    Zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
    KalmanResid(arma::span(),arma::span(i-1,i-1)) = EpsilonYBetaPlus(arma::span(),arma::span(i-1,i-1)) - Zt*BetaFiltered(arma::span(),arma::span(i-1,i-1));
    StateCovPredicted = Zt*StateCovFiltered*trans(Zt) + SigmaDraw;
    invStateCovPredicted = inv_sympd(StateCovPredicted);
    FStore(arma::span((i-1)*M,(i*M)-1),arma::span()) = invStateCovPredicted;
    KalmanGain = StateCovFiltered*trans(Zt)*invStateCovPredicted;
    Lt = arma::eye(K,K) - KalmanGain*Zt;
    LStore(arma::span((i-1)*K,(i*K)-1),arma::span()) = Lt;
    BetaFiltered(arma::span(),arma::span(i,i)) = BetaFiltered(arma::span(),arma::span(i-1,i-1)) + KalmanGain*KalmanResid(arma::span(),arma::span(i-1,i-1));
    StateCovFiltered = StateCovFiltered*trans(Lt) + QDraw;
  }
  //
  for(j=kT; j>=1; j--){
    Zt = Z(arma::span((j-1)*M,(j*M)-1),arma::span());
    Lt = LStore(arma::span((j-1)*K,(j*K)-1),arma::span());
    invStateCovPredicted = FStore(arma::span((j-1)*M,(j*M)-1),arma::span());
    rt(arma::span(),arma::span(j-1,j-1)) = trans(Zt)*invStateCovPredicted*KalmanResid(arma::span(),arma::span(j-1,j-1)) + trans(Lt)*rt(arma::span(),arma::span(j,j));
  }
  //
  for(i=1; i<=kT; i++) {
    BetaSmoothed(arma::span(),arma::span(i,i)) = BetaSmoothed(arma::span(),arma::span(i-1,i-1)) + QDraw*rt(arma::span(),arma::span(i,i));
  }
  //
  BetaTilde = BetaSmoothed + BetaPlus;
  //
  for(i=1; i <=(kT+1); i++){
    BtDraw(arma::span(),arma::span(i-1,i-1)) = BetaTilde(arma::span(),arma::span(i-1,i-1)) + B0Draw;
  }
  arma::mat DeltaBeta = trans(BtDraw(arma::span(),arma::span(1,kT))) - trans(BtDraw(arma::span(),arma::span(0,kT-1)));
  arma::mat SSEQ(K,K);
  SSEQ.zeros();
  for(i=1; i<=kT; i++){
    SSEQ += trans(DeltaBeta(arma::span(i-1,i-1),arma::span()))*DeltaBeta(arma::span(i-1,i-1),arma::span());
  }
  //
  arma::mat Qinv = inv_sympd(SSEQ + QPr);
  arma::mat Krand1 = arma::randn(Qinv.n_rows,kT+QVPr);
  arma::mat IWA = trans(arma::chol(Qinv))*Krand1;
  arma::mat QinvDraw = IWA*trans(IWA);
  QDraw = inv_sympd(QinvDraw);
  QChol = arma::chol(QDraw);
  //
  Epsilon.zeros();
  for(i=1; i<=kT; i++){
    Epsilon(arma::span(),arma::span(i-1,i-1)) = y(arma::span(),arma::span(i-1,i-1)) - Z(arma::span((i-1)*M,(i*M)-1),arma::span())*BtDraw(arma::span(),arma::span(i-1,i-1));
  }
  arma::mat SSES(M,M);
  SSES.zeros();
  for(i=1; i<=kT; i++){
    SSES += Epsilon(arma::span(),arma::span(i-1,i-1))*trans(Epsilon(arma::span(),arma::span(i-1,i-1)));
  }
  //
  arma::mat Sigmainv = inv_sympd(SSES + SPr);
  arma::mat Krand2 = arma::randn(Sigmainv.n_rows,kT+SVPr);
  arma::mat IWA2 = trans(arma::chol(Sigmainv))*Krand2;
  invSDraw = IWA2*trans(IWA2);
  SigmaDraw = inv_sympd(invSDraw);
  SChol = arma::chol(SigmaDraw);
  //
  //
  //
  for(kreps=1; kreps<=(keep+burnin); kreps++){
    //
    BVPt.zeros();
    BetaNumerator.zeros();
    Zt.zeros();
    //
    for(i=1; i<=kT; i++){
      Zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
      BVPt += trans(Zt)*invSDraw*Zt;
      BetaNumerator += trans(Zt)*invSDraw*(y(arma::span(),arma::span(i-1,i-1)) - Zt*BetaTilde(arma::span(),arma::span(i-1,i-1)));
    }
    BVPt = inv_sympd(BVPt + invB0VPr);
    KBRand = arma::randn(K,1);
    B0Draw = BVPt*(invB0VPr*B0Pr + BetaNumerator) + trans(arma::chol(BVPt))*KBRand;
    //
    Epsilon.zeros();
    for(i=1; i<=kT; i++){
      Epsilon(arma::span(),arma::span(i-1,i-1)) = y(arma::span(),arma::span(i-1,i-1)) - Z(arma::span((i-1)*M,(i*M)-1),arma::span())*B0Draw;
    }
    //
    ResidPlus.zeros();
    KWRand1 = arma::randn(M,1);
    KWRand2 = arma::randn(K,1);
    //
    for(i=1; i<=kT; i++){
      ResidPlus(arma::span((i-1)*(M+K),(i-1)*(M+K)+M-1),arma::span(0,0)) = trans(SChol)*KWRand1;
      ResidPlus(arma::span((i-1)*(M+K)+M,(i-1)*(M+K)+(M+K)-1),arma::span(0,0)) = trans(QChol)*KWRand2;
      //
      KWRand1 = arma::randn(M,1);
      KWRand2 = arma::randn(K,1);
    }
    //
    Zt.zeros();
    BetaPlus.zeros();
    YBetaPlus.zeros();
    //
    for(i=1; i<=kT; i++) {
      Zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
      YBetaPlus(arma::span(),arma::span(i-1,i-1)) = Zt*BetaPlus(arma::span(),arma::span(i-1,i-1)) + ResidPlus(arma::span((i-1)*(M+K),(i-1)*(M+K)+M-1),arma::span());
      BetaPlus(arma::span(),arma::span(i,i)) = BetaPlus(arma::span(),arma::span(i-1,i-1)) + ResidPlus(arma::span((i-1)*(M+K) + M,(i*(M+K))-1),arma::span());
    }
    //
    EpsilonYBetaPlus.zeros();
    EpsilonYBetaPlus = Epsilon - YBetaPlus;
    //
    LStore.zeros();
    FStore.zeros();
    BetaFiltered.zeros();
    KalmanResid.zeros();
    StateCovFiltered.zeros();
    rt.zeros();
    BetaSmoothed.zeros();
    //
    Zt = Z(arma::span(0,M-1),arma::span());
    KalmanResid(arma::span(),arma::span(0,0)) = EpsilonYBetaPlus(arma::span(),arma::span(0,0)) - Zt*BetaFiltered(arma::span(),arma::span(0,0));
    StateCovPredicted = Zt*StateCovFiltered*trans(Zt) + SigmaDraw;
    invStateCovPredicted = inv_sympd(StateCovPredicted);
    FStore(arma::span(0,M-1),arma::span()) = invStateCovPredicted;
    KalmanGain = StateCovFiltered*trans(Zt)*invStateCovPredicted;
    Lt = arma::eye(K,K) - KalmanGain*Zt;
    LStore(arma::span(0,K-1),arma::span()) = Lt;
    BetaFiltered(arma::span(),arma::span(1,1)) = BetaFiltered(arma::span(),arma::span(0,0)) + KalmanGain*KalmanResid(arma::span(),arma::span(0,0));
    StateCovFiltered = StateCovFiltered*trans(Lt) + QDraw;
    //
    for(i=2; i<=kT; i++) {
      Zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
      KalmanResid(arma::span(),arma::span(i-1,i-1)) = EpsilonYBetaPlus(arma::span(),arma::span(i-1,i-1)) - Zt*BetaFiltered(arma::span(),arma::span(i-1,i-1));
      StateCovPredicted = Zt*StateCovFiltered*trans(Zt) + SigmaDraw;
      invStateCovPredicted = inv_sympd(StateCovPredicted);
      FStore(arma::span((i-1)*M,(i*M)-1),arma::span()) = invStateCovPredicted;
      KalmanGain = StateCovFiltered*trans(Zt)*invStateCovPredicted;
      Lt = arma::eye(K,K) - KalmanGain*Zt;
      LStore(arma::span((i-1)*K,(i*K)-1),arma::span()) = Lt;
      BetaFiltered(arma::span(),arma::span(i,i)) = BetaFiltered(arma::span(),arma::span(i-1,i-1)) + KalmanGain*KalmanResid(arma::span(),arma::span(i-1,i-1));
      StateCovFiltered = StateCovFiltered*trans(Lt) + QDraw;
    }
    //
    for(j=kT; j>=1; j--){
      Zt = Z(arma::span((j-1)*M,(j*M)-1),arma::span());
      Lt = LStore(arma::span((j-1)*K,(j*K)-1),arma::span());
      invStateCovPredicted = FStore(arma::span((j-1)*M,(j*M)-1),arma::span());
      rt(arma::span(),arma::span(j-1,j-1)) = trans(Zt)*invStateCovPredicted*KalmanResid(arma::span(),arma::span(j-1,j-1)) + trans(Lt)*rt(arma::span(),arma::span(j,j));
    }
    //
    for(i=1; i<=kT; i++) {
      BetaSmoothed(arma::span(),arma::span(i,i)) = BetaSmoothed(arma::span(),arma::span(i-1,i-1)) + QDraw*rt(arma::span(),arma::span(i,i));
    }
    //
    BetaTilde = BetaSmoothed + BetaPlus;
    //
    for(i=1; i <=(kT+1); i++){
      BtDraw(arma::span(),arma::span(i-1,i-1)) = BetaTilde(arma::span(),arma::span(i-1,i-1)) + B0Draw;
    }
    DeltaBeta = trans(BtDraw(arma::span(),arma::span(1,kT))) - trans(BtDraw(arma::span(),arma::span(0,kT-1)));
    SSEQ.zeros();
    for(i=1; i<=kT; i++){
      SSEQ += trans(DeltaBeta(arma::span(i-1,i-1),arma::span()))*DeltaBeta(arma::span(i-1,i-1),arma::span());
    }
    //
    Qinv = inv_sympd(SSEQ + QPr);
    Krand1 = arma::randn(Qinv.n_rows,kT+QVPr);
    IWA = trans(arma::chol(Qinv))*Krand1;
    QinvDraw = IWA*trans(IWA);
    QDraw = inv_sympd(QinvDraw);
    QChol = arma::chol(QDraw);
    //
    Epsilon.zeros();
    for(i=1; i<=kT; i++){
      Epsilon(arma::span(),arma::span(i-1,i-1)) = y(arma::span(),arma::span(i-1,i-1)) - Z(arma::span((i-1)*M,(i*M)-1),arma::span())*BtDraw(arma::span(),arma::span(i-1,i-1));
    }
    SSES.zeros();
    for(i=1; i<=kT; i++){
      SSES += Epsilon(arma::span(),arma::span(i-1,i-1))*trans(Epsilon(arma::span(),arma::span(i-1,i-1)));
    }
    //
    Sigmainv = inv_sympd(SSES + SPr);
    Krand2 = arma::randn(Sigmainv.n_rows,kT+SVPr);
    IWA2 = trans(arma::chol(Sigmainv))*Krand2;
    invSDraw = IWA2*trans(IWA2);
    SigmaDraw = inv_sympd(invSDraw);
    SChol = arma::chol(SigmaDraw);
    //
    if(kreps > burnin){
      BetaDraws(arma::span(),arma::span(),arma::span(kreps-burnin-1,kreps-burnin-1)) = BtDraw;
      QDraws(arma::span(),arma::span(),arma::span(kreps-burnin-1,kreps-burnin-1)) = QDraw;
      SDraws(arma::span(),arma::span(),arma::span(kreps-burnin-1,kreps-burnin-1)) = SigmaDraw;
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("BetaDraws") = BetaDraws,Rcpp::Named("QDraws") = QDraws,Rcpp::Named("SDraws") = SDraws);
}

SEXP BVARTVPIRFs( SEXP mM, SEXP mK, SEXP mkreps, SEXP mirfperiods, 
                  SEXP mBDs, SEXP mSDs )
{
  int M = as<int>(mM);
  int K = as<int>(mK);
  int kreps = as<int>(mkreps);
  int irfperiods = as<int>(mirfperiods);
  //
  int kcons = 1;
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
