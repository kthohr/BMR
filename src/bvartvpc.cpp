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

#include "bvartvpc.h"
using namespace Rcpp;

SEXP tsprior( SEXP mY, SEXP mtau, SEXP mM, SEXP mK, SEXP mp )
{
    try {
        arma::mat Y = as<arma::mat>(mY);
        int tau = as<int>(mtau);
        int M = as<int>(mM);
        int K = as<int>(mK);
        int p = as<int>(mp);
        //
        arma::mat Z = arma::zeros<arma::mat>(M*tau,K);
        //
        arma::mat BetaDenominator = arma::zeros<arma::mat>(K,K);
        arma::mat BetaNumerator = arma::zeros<arma::mat>(K,1);
        //
        int i,j,k;
        //
        //
        //
        arma::mat Ytau = arma::trans(Y.rows(p,tau+p-1));
        //
        for(i=(p+1); i<=(tau+p); i++){
            k = i - p - 1;
            Z(arma::span(k*M,(k+1)*M-1),arma::span(0,M-1)) = arma::eye(M,M);
            for(j=1; j<=p; j++){
                Z(arma::span(k*M, (k+1)*M-1),arma::span(M+((M*M)*(j-1)), M-1+((M*M)*j))) = arma::kron(arma::eye(M,M), Y.row(i-1-j));
            }
        }
        //
        arma::mat zt = arma::zeros<arma::mat>(M,K);
        arma::mat Beta = arma::zeros<arma::mat>(K,1);
        //
        for(i=1;i<=tau;i++){
            zt = Z.rows((i-1)*M,(i*M)-1);
            BetaDenominator += zt.t()*zt;
            BetaNumerator += zt.t()*Ytau.col(i-1);
        }
        Beta = arma::inv(BetaDenominator)*BetaNumerator;
        //
        arma::mat SSE = arma::zeros<arma::mat>(M,M);
        //
        for(i=1;i<=tau;i++){
            zt = Z(arma::span((i-1)*M,(i*M)-1),arma::span());
            SSE += (Ytau.col(i-1) - zt*Beta)*arma::trans(Ytau.col(i-1) - zt*Beta);
        }
        //
        double tau2 = tau - (M*p + 1);
        arma::mat SPrOLS = SSE*(1/tau2);
        arma::mat SPrOLSinv = arma::inv(SPrOLS);
        //
        arma::mat BVPrOLS = arma::zeros<arma::mat>(K,K);
        for(i=1;i<=tau;i++){
            zt = Z.rows((i-1)*M,(i*M)-1);
            BVPrOLS += zt.t()*SPrOLSinv*zt;
        }
        BVPrOLS = arma::inv(BVPrOLS);
        //
        return Rcpp::List::create(Rcpp::Named("BetaOLS") = Beta,Rcpp::Named("BVPrOLS") = BVPrOLS, Rcpp::Named("SigmaOLS") = SPrOLS);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARTVP Prior C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP BVARTVPReps( SEXP my, SEXP mZ, SEXP mM, SEXP mK, SEXP mkT, SEXP mkeep, 
                  SEXP mburnin, SEXP mB0Pr, SEXP mB0VPr, SEXP minvB0VPr, SEXP mQPr,
                  SEXP mQVPr, SEXP mSPr, SEXP mSVPr, SEXP mQDraw, SEXP mQChol,
                  SEXP mSDraw, SEXP minvSDraw )
{
    try {
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
        arma::mat BtDraw = arma::zeros<arma::mat>(K,(kT+1));
        arma::mat BetaTilde = arma::zeros<arma::mat>(K,(kT+1));
        //
        arma::mat SChol = arma::chol(SigmaDraw);
        //
        arma::mat BVPt = arma::zeros<arma::mat>(K,K);
        arma::mat BetaNumerator = arma::zeros<arma::mat>(K,1);
        arma::mat Zt = arma::zeros<arma::mat>(M,K);
        //
        for(i=1; i<=kT; i++){
            Zt = Z.rows((i-1)*M,(i*M)-1);
            BVPt += Zt.t()*invSDraw*Zt;
            BetaNumerator += Zt.t()*invSDraw*(y.col(i-1) - Zt*BetaTilde.col(i-1));
        }
        BVPt = arma::inv_sympd(BVPt + invB0VPr);
        arma::mat KBRand = arma::randn(K,1);
        arma::mat B0Draw = BVPt*(invB0VPr*B0Pr + BetaNumerator) + arma::trans(arma::chol(BVPt))*KBRand;
        //
        arma::mat Epsilon = arma::zeros<arma::mat>(M,kT);
        arma::mat BetaPlus = arma::zeros<arma::mat>(K,kT+1);
        arma::mat YBetaPlus = arma::zeros<arma::mat>(M,kT);
        //
        arma::mat KWRand1;
        arma::mat KWRand2;
        //
        for(i=1; i<=kT; i++){
            Epsilon.col(i-1) = y.col(i-1) - Z.rows((i-1)*M,(i*M)-1)*B0Draw;
            //
            KWRand1 = arma::randn(M,1);
            KWRand2 = arma::randn(K,1);
            //
            Zt = Z.rows((i-1)*M,(i*M)-1);
            //
            YBetaPlus.col(i-1) = Zt*BetaPlus.col(i-1) + SChol.t()*KWRand1;
            BetaPlus.col(i) = BetaPlus.col(i-1) + QChol.t()*KWRand2;
        }
        //
        arma::mat EpsilonYBetaPlus = Epsilon - YBetaPlus;
        //
        arma::mat LStore = arma::zeros<arma::mat>(K*kT,K);
        arma::mat FStore = arma::zeros<arma::mat>(M*kT,M);
        arma::mat BetaFiltered = arma::zeros<arma::mat>(K,kT+1);
        arma::mat KalmanResid = arma::zeros<arma::mat>(M,kT);
        arma::mat StateCovFiltered = arma::zeros<arma::mat>(K,K);
        //
        arma::mat rt = arma::zeros<arma::mat>(K,kT+1);
        arma::mat BetaSmoothed = arma::zeros<arma::mat>(K,kT+1);
        //
        Zt.zeros();
        //
        Zt = Z.rows(0,M-1);
        KalmanResid.col(0) = EpsilonYBetaPlus.col(0) - Zt*BetaFiltered.col(0);
        arma::mat StateCovPredicted = Zt*StateCovFiltered*Zt.t() + SigmaDraw;
        arma::mat invStateCovPredicted = arma::inv_sympd(StateCovPredicted);
        FStore.rows(0,M-1) = invStateCovPredicted;
        arma::mat KalmanGain = StateCovFiltered*Zt.t()*invStateCovPredicted;
        arma::mat Lt = arma::eye(K,K) - KalmanGain*Zt;
        LStore.rows(0,K-1) = Lt;
        BetaFiltered.col(1) = BetaFiltered.col(0) + KalmanGain*KalmanResid.col(0);
        StateCovFiltered = StateCovFiltered*Lt.t() + QDraw;
        //
        for(i=2; i<=kT; i++){
            Zt = Z.rows((i-1)*M,(i*M)-1);
            KalmanResid.col(i-1) = EpsilonYBetaPlus.col(i-1) - Zt*BetaFiltered.col(i-1);
            StateCovPredicted = Zt*StateCovFiltered*Zt.t() + SigmaDraw;
            invStateCovPredicted = arma::inv_sympd(StateCovPredicted);
            FStore.rows((i-1)*M,(i*M)-1) = invStateCovPredicted;
            KalmanGain = StateCovFiltered*Zt.t()*invStateCovPredicted;
            Lt = arma::eye(K,K) - KalmanGain*Zt;
            LStore.rows((i-1)*K,(i*K)-1) = Lt;
            BetaFiltered.col(i) = BetaFiltered.col(i-1) + KalmanGain*KalmanResid.col(i-1);
            StateCovFiltered = StateCovFiltered*Lt.t() + QDraw;
        }
        //
        for(j=kT; j>=1; j--){
            Zt = Z.rows((j-1)*M,(j*M)-1);
            Lt = LStore.rows((j-1)*K,(j*K)-1);
            invStateCovPredicted = FStore.rows((j-1)*M,(j*M)-1);
            rt.col(j-1) = Zt.t()*invStateCovPredicted*KalmanResid.col(j-1) + Lt.t()*rt.col(j);
        }
        //
        for(i=1; i<=kT; i++){
            BetaSmoothed.col(i) = BetaSmoothed.col(i-1) + QDraw*rt.col(i);
        }
        //
        BetaTilde = BetaSmoothed + BetaPlus;
        //
        BtDraw = BetaTilde + arma::repmat(B0Draw,1,kT+1);
        //
        arma::mat DeltaBeta = arma::trans(BtDraw.cols(1,kT)) - arma::trans(BtDraw.cols(0,kT-1));
        arma::mat SSEQ = DeltaBeta.t()*DeltaBeta;
        //
        arma::mat Qinv = arma::inv_sympd(SSEQ + QPr);
        arma::mat Krand1 = arma::randn(Qinv.n_rows,kT+QVPr);
        arma::mat IWA = arma::trans(arma::chol(Qinv))*Krand1;
        arma::mat QinvDraw = IWA*arma::trans(IWA);
        QDraw = arma::inv_sympd(QinvDraw);
        QChol = arma::chol(QDraw);
        //
        Epsilon.zeros();
        for(i=1; i<=kT; i++){
            Epsilon.col(i-1) = y.col(i-1) - Z.rows((i-1)*M,(i*M)-1)*BtDraw.col(i-1);
        }
        arma::mat SSES = Epsilon*Epsilon.t();
        //
        arma::mat Sigmainv = arma::inv_sympd(SSES + SPr);
        arma::mat Krand2 = arma::randn(Sigmainv.n_rows,kT+SVPr);
        arma::mat IWA2 = arma::trans(arma::chol(Sigmainv))*Krand2;
        invSDraw = IWA2*trans(IWA2);
        SigmaDraw = arma::inv_sympd(invSDraw);
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
                Zt = Z.rows((i-1)*M,(i*M)-1);
                BVPt += Zt.t()*invSDraw*Zt;
                BetaNumerator += Zt.t()*invSDraw*(y.col(i-1) - Zt*BetaTilde.col(i-1));
            }
            BVPt = arma::inv_sympd(BVPt + invB0VPr);
            KBRand = arma::randn(K,1);
            B0Draw = BVPt*(invB0VPr*B0Pr + BetaNumerator) + arma::trans(arma::chol(BVPt))*KBRand;
            //
            Epsilon.zeros();
            Zt.zeros();
            BetaPlus.zeros();
            YBetaPlus.zeros();
            //
            for(i=1; i<=kT; i++){
                Epsilon.col(i-1) = y.col(i-1) - Z.rows((i-1)*M,(i*M)-1)*B0Draw;
                //
                KWRand1 = arma::randn(M,1);
                KWRand2 = arma::randn(K,1);
                //
                Zt = Z.rows((i-1)*M,(i*M)-1);
                //
                YBetaPlus.col(i-1) = Zt*BetaPlus.col(i-1) + SChol.t()*KWRand1;
                BetaPlus.col(i) = BetaPlus.col(i-1) + QChol.t()*KWRand2;
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
            Zt = Z.rows(0,M-1);
            KalmanResid.col(0) = EpsilonYBetaPlus.col(0) - Zt*BetaFiltered.col(0);
            StateCovPredicted = Zt*StateCovFiltered*Zt.t() + SigmaDraw;
            invStateCovPredicted = arma::inv_sympd(StateCovPredicted);
            FStore.rows(0,M-1) = invStateCovPredicted;
            KalmanGain = StateCovFiltered*Zt.t()*invStateCovPredicted;
            Lt = arma::eye(K,K) - KalmanGain*Zt;
            LStore.rows(0,K-1) = Lt;
            BetaFiltered.col(1) = BetaFiltered.col(0) + KalmanGain*KalmanResid.col(0);
            StateCovFiltered = StateCovFiltered*Lt.t() + QDraw;
            //
            for(i=2; i<=kT; i++) {
                Zt = Z.rows((i-1)*M,(i*M)-1);
                KalmanResid.col(i-1) = EpsilonYBetaPlus.col(i-1) - Zt*BetaFiltered.col(i-1);
                StateCovPredicted = Zt*StateCovFiltered*Zt.t() + SigmaDraw;
                invStateCovPredicted = arma::inv_sympd(StateCovPredicted);
                FStore.rows((i-1)*M,(i*M)-1) = invStateCovPredicted;
                KalmanGain = StateCovFiltered*Zt.t()*invStateCovPredicted;
                Lt = arma::eye(K,K) - KalmanGain*Zt;
                LStore.rows((i-1)*K,(i*K)-1) = Lt;
                BetaFiltered.col(i) = BetaFiltered.col(i-1) + KalmanGain*KalmanResid.col(i-1);
                StateCovFiltered = StateCovFiltered*Lt.t() + QDraw;
            }
            //
            for(j=kT; j>=1; j--){
                Zt = Z.rows((j-1)*M,(j*M)-1);
                Lt = LStore.rows((j-1)*K,(j*K)-1);
                invStateCovPredicted = FStore.rows((j-1)*M,(j*M)-1);
                rt.col(j-1) = Zt.t()*invStateCovPredicted*KalmanResid.col(j-1) + Lt.t()*rt.col(j);
            }
            //
            for(i=1; i<=kT; i++) {
                BetaSmoothed.col(i) = BetaSmoothed.col(i-1) + QDraw*rt.col(i);
            }
            //
            BetaTilde = BetaSmoothed + BetaPlus;
            //
            BtDraw = BetaTilde + arma::repmat(B0Draw,1,kT+1);
            //
            DeltaBeta = arma::trans(BtDraw.cols(1,kT)) - arma::trans(BtDraw.cols(0,kT-1));
            arma::mat SSEQ = DeltaBeta.t()*DeltaBeta;
            //
            Qinv = arma::inv_sympd(SSEQ + QPr);
            Krand1 = arma::randn(Qinv.n_rows,kT+QVPr);
            IWA = arma::trans(arma::chol(Qinv))*Krand1;
            QinvDraw = IWA*IWA.t();
            QDraw = arma::inv_sympd(QinvDraw);
            QChol = arma::chol(QDraw);
            //
            Epsilon.zeros();
            for(i=1; i<=kT; i++){
                Epsilon.col(i-1) = y.col(i-1) - Z.rows((i-1)*M,(i*M)-1)*BtDraw.col(i-1);
            }
            SSES = Epsilon*Epsilon.t();
            //
            Sigmainv = arma::inv_sympd(SSES + SPr);
            Krand2 = arma::randn(Sigmainv.n_rows,kT+SVPr);
            IWA2 = arma::trans(arma::chol(Sigmainv))*Krand2;
            invSDraw = IWA2*trans(IWA2);
            SigmaDraw = arma::inv_sympd(invSDraw);
            SChol = arma::chol(SigmaDraw);
            //
            if(kreps > burnin){
                BetaDraws.slice(kreps-burnin-1) = BtDraw;
                QDraws.slice(kreps-burnin-1) = QDraw;
                SDraws.slice(kreps-burnin-1) = SigmaDraw;
            }
        }
        //
        return Rcpp::List::create(Rcpp::Named("BetaDraws") = BetaDraws,Rcpp::Named("QDraws") = QDraws,Rcpp::Named("SDraws") = SDraws);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARTVP Gibbs C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP BVARTVPIRFs( SEXP mM, SEXP mK, SEXP mkreps, SEXP mirfperiods, 
                  SEXP mBDs, SEXP mSDs )
{
    try {
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
        arma::mat BetaT = arma::zeros<arma::mat>(K-kcons,M);
        arma::mat shock = arma::zeros<arma::mat>(M,M);
        arma::mat ShockO = arma::zeros<arma::mat>(K-M-kcons,M);
        arma::mat ShockB = arma::zeros<arma::mat>(K-kcons,M);
        arma::mat Pow = arma::zeros<arma::mat>(M,M);
        //
        for (j=1;j<=kreps;j++) {
            BetaT = BetaDraws(arma::span(kcons,K-1),arma::span(),arma::span(j-1,j-1));
            BetaT = arma::trans(BetaT);
            //
            shock = SigmaDraws.slice(j-1);
            shock = arma::trans(arma::chol(shock));
            //
            ShockB.zeros();
            ShockB.rows(0,M-1) = shock;
            //
            ImpStore.slice((j-1)*irfperiods) = shock;
            //
            for(i=2; i<=irfperiods; i++){
                Pow = BetaT*ShockB;
                ImpStore.slice((j-1)*irfperiods + (i-1)) = Pow;
                if(K > M+kcons){
                    ShockB.rows(M,K-1-kcons) = ShockB.rows(0,K-M-1-kcons);
                }
                ShockB.rows(0,M-1) = Pow;
            }
        }
        //
        return Rcpp::List::create(Rcpp::Named("ImpStore") = ImpStore);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARTVP IRF C++ exception (unknown reason)" );
    }
    return R_NilValue;
}