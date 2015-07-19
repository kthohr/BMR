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

#include "bvarsc.h"
using namespace Rcpp;

SEXP SBVARReps(SEXP mX, SEXP mY, SEXP md, SEXP mdX, SEXP myd, SEXP mZd, SEXP mPsiPr,
               SEXP minvPsiVPr, SEXP mBPr, SEXP mBeta, SEXP minvBVPr, SEXP mSigma,
               SEXP mSigmaML, SEXP mgamma, SEXP mTp, SEXP mM, SEXP mp, SEXP mburnin,
               SEXP mkreps)
{
    try {
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
        //
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
        arma::mat D = arma::zeros<arma::mat>(Tp,p+1);
        D.col(0) = d;
        D.cols(1,p) = dX;
        //
        arma::mat UFixed = arma::zeros<arma::mat>(M + p*M*q, M*q);
        UFixed.rows(0,M-1) = arma::eye(M,M);
        arma::mat U = UFixed;
        //
        // Psi
        //
        arma::mat xi = Y - X*BetaKT;
        U = UFixed;
        arma::mat BetaT = BetaKT.t();
        for(j=1;j<=p;j++){
            U.rows(M*q*j,M*q*(j+1)-1) = BetaT.cols(M*(j-1),M*j-1);
        }
        //
        arma::mat invSigma = arma::inv_sympd(Sigma);
        arma::mat invPsiVPt = U.t()*arma::kron(D.t()*D,invSigma)*U + invPsiVPr;
        arma::colvec StackedZD = arma::vectorise(invSigma*xi.t()*D);
        arma::mat PsiVPt = arma::inv_sympd(invPsiVPt);
        arma::mat vecPsiPt = PsiVPt*(U.t()*StackedZD + invPsiVPr*PsiPr.t());
        KPsiRand = arma::randn(M*q);
        arma::mat vecPsi = vecPsiPt + arma::trans(arma::chol(PsiVPt))*KPsiRand;
        arma::mat Psi = arma::reshape(vecPsi,1,M);
        //
        // Beta
        yd = Y - d*Psi;
        Zd = X + dX*arma::kron(arma::eye(p,p),Psi);
        arma::mat BetaVPt = arma::inv_sympd(Zd.t()*Zd + invBVPr);
        arma::mat BetaPt = BetaVPt*(Zd.t()*yd + invBVPr*BPr);
        KBetaRand = arma::randn(M*M*p);
        arma::colvec vecBetaT = arma::vectorise(BetaPt);
        arma::mat vecBeta = vecBetaT + arma::trans(arma::chol(arma::kron(Sigma,BetaVPt)))*KBetaRand;
        arma::mat Beta = arma::reshape(vecBeta,M*p,M);
        //
        // Sigma
        arma::mat Epsilon = yd - Zd*Beta;
        arma::mat B = arma::trans(Beta - BPr)*invBVPr*(Beta - BPr);
        Sigma = Epsilon.t()*Epsilon + SigmaML + B;
        KSigmaRand = arma::randn(Sigma.n_rows,Tp + p*M + gamma);
        Sigma = arma::trans(arma::chol(arma::inv_sympd(Sigma)))*KSigmaRand;
        Sigma = arma::inv_sympd(Sigma*Sigma.t());
        //
        //
        //
        for(i=1; i<=burnin; i++){
            xi = Y - X*Beta;
            U = UFixed;
            arma::mat BetaT = Beta.t();
            //
            for(j=1;j<=p;j++){
                U.rows(M*q*j, M*q*(j+1)-1) = BetaT.cols(M*(j-1), M*j-1);
            }
            //
            invSigma = arma::inv_sympd(Sigma);
            invPsiVPt = U.t()*arma::kron(D.t()*D,invSigma)*U + invPsiVPr;
            StackedZD = arma::vectorise(invSigma*xi.t()*D);;
            PsiVPt = arma::inv_sympd(invPsiVPt);
            vecPsiPt = PsiVPt*(U.t()*StackedZD + invPsiVPr*PsiPr.t());
            KPsiRand = arma::randn(M*q);
            vecPsi = vecPsiPt + arma::trans(arma::chol(PsiVPt))*KPsiRand;
            Psi = arma::reshape(vecPsi,1,M);
            //
            // Beta
            yd = Y - d*Psi;
            Zd = X + dX*arma::kron(arma::eye(p,p),Psi);
            BetaVPt = arma::inv_sympd(Zd.t()*Zd + invBVPr);
            BetaPt = BetaVPt*(Zd.t()*yd + invBVPr*BPr);
            KBetaRand = arma::randn(M*M*p);
            vecBetaT = arma::vectorise(BetaPt);
            vecBeta = vecBetaT + arma::trans(arma::chol(arma::kron(Sigma,BetaVPt)))*KBetaRand;
            Beta = arma::reshape(vecBeta,M*p,M);
            //
            // Sigma
            Epsilon = yd - Zd*Beta;
            B = arma::trans(Beta - BPr)*invBVPr*(Beta - BPr);
            Sigma = Epsilon.t()*Epsilon + SigmaML + B;
            KSigmaRand = arma::randn(Sigma.n_rows,Tp + p*M + gamma);
            Sigma = arma::trans(arma::chol(arma::inv_sympd(Sigma)))*KSigmaRand;
            Sigma = arma::inv_sympd(Sigma*Sigma.t());
        }
        //
        for(i=1; i<=kreps; i++){
            xi = Y - X*Beta;
            U = UFixed;
            arma::mat BetaT = Beta.t();
            for(j=1; j<=p; j++){
                U.rows(M*q*j, M*q*(j+1)-1) = BetaT.cols(M*(j-1), M*j-1);
            }
            //
            invSigma = arma::inv_sympd(Sigma);
            invPsiVPt = U.t()*arma::kron(D.t()*D,invSigma)*U + invPsiVPr;
            StackedZD = arma::vectorise(invSigma*xi.t()*D);;
            PsiVPt = arma::inv_sympd(invPsiVPt);
            vecPsiPt = PsiVPt*(U.t()*StackedZD + invPsiVPr*PsiPr.t());
            KPsiRand = arma::randn(M*q);
            vecPsi = vecPsiPt + arma::trans(arma::chol(PsiVPt))*KPsiRand;
            Psi = arma::reshape(vecPsi,1,M);
            //
            // Beta
            yd = Y - d*Psi;
            Zd = X + dX*arma::kron(arma::eye(p,p),Psi);
            BetaVPt = arma::inv_sympd(Zd.t()*Zd + invBVPr);
            BetaPt = BetaVPt*(Zd.t()*yd + invBVPr*BPr);
            KBetaRand = arma::randn(M*M*p);
            vecBetaT = arma::vectorise(BetaPt);
            vecBeta = vecBetaT + arma::trans(arma::chol(arma::kron(Sigma,BetaVPt)))*KBetaRand;
            Beta = arma::reshape(vecBeta,M*p,M);
            //
            // Sigma
            Epsilon = yd - Zd*Beta;
            B = arma::trans(Beta - BPr)*invBVPr*(Beta - BPr);
            Sigma = Epsilon.t()*Epsilon + SigmaML + B;
            KSigmaRand = arma::randn(Sigma.n_rows,Tp + p*M + gamma);
            Sigma = arma::trans(arma::chol(arma::inv_sympd(Sigma)))*KSigmaRand;
            Sigma = arma::inv_sympd(Sigma*Sigma.t());
            //
            PsiDraws.slice(i-1) = Psi;
            BetaDraws.slice(i-1) = Beta;
            SigmaDraws.slice(i-1) = Sigma;
        }
        //
        return Rcpp::List::create(Rcpp::Named("Psi") = PsiDraws,Rcpp::Named("Beta") = BetaDraws,Rcpp::Named("Sigma") = SigmaDraws);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARS Gibbs C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP SBVARIRFs( SEXP mM, SEXP mK, SEXP mkreps, SEXP mirfperiods, SEXP mBDs, SEXP mSDs )
{
    try {
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
        arma::mat BetaT = arma::zeros<arma::mat>(K,M);
        arma::mat shock = arma::zeros<arma::mat>(M,M);
        arma::mat ShockO = arma::zeros<arma::mat>(K-M,M);
        arma::mat ShockB = arma::zeros<arma::mat>(K,M);
        arma::mat Pow = arma::zeros<arma::mat>(M,M);
        //
        for(j=1; j<=kreps; j++) {
            BetaT = BetaDraws.slice(j-1);
            BetaT = BetaT.t();
            shock = SigmaDraws.slice(j-1);
            shock = arma::trans(arma::chol(shock));
            ShockB.zeros();
            ShockB.rows(0,M-1) = shock;
            ImpStore.slice((j-1)*irfperiods) = shock;
            for(i=2; i<=irfperiods; i++){
                Pow = BetaT*ShockB;
                ImpStore.slice((j-1)*irfperiods + (i-1)) = Pow;
                if (K > M){
                    ShockB.rows(M,K-1) = ShockB.rows(0,K-M-1);
                }
                ShockB.rows(0,M-1) = Pow;
            }
        }
        //
        return Rcpp::List::create(Rcpp::Named("ImpStore") = ImpStore);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARS IRF C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP bvarsforecast( SEXP mY0, SEXP mdX, SEXP mM, SEXP mp, SEXP mK, SEXP mruns, SEXP mperiods,
                    SEXP minclshocks, SEXP mPDs, SEXP mBDs, SEXP mSDs )
{
    try {
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
        arma::mat PT = PsiDraws.slice(0);
        arma::mat BT = BetaDraws.slice(0);
        arma::mat ST = SigmaDraws.slice(0);
        //
        arma::mat Shock = arma::trans(arma::chol(ST));
        //
        arma::mat fY = arma::zeros<arma::mat>(periods,M);
        arma::mat kY = Y0;
        arma::vec Krand = arma::randn(M);
        arma::mat jshock = arma::zeros<arma::mat>(1,M);
        //
        int i,j;
        //
        if(inclshocks > 0){
            for(i=1; i<=runs; i++){
                PT = PsiDraws.slice(i-1);
                BT = BetaDraws.slice(i-1);
                ST = SigmaDraws.slice(i-1);
                //
                Shock = arma::trans(arma::chol(ST));
                //
                fY.zeros();
                kY = Y0;
                //
                for(j=1; j<=periods; j++){
                    Krand = arma::randn(M);
                    jshock = arma::trans(Shock*Krand);
                    //
                    fY.row(j-1) = kY*BT + PT - dX*(arma::kron(arma::eye(p,p),PT))*BT + jshock;
                    if (K > M){
                        kY(0,arma::span(M,K-1)) = kY(0,arma::span(0,K-M-1));
                    }
                    kY(0,arma::span(0,M-1)) = fY.row(j-1);
                }
                Forecasts.slice(i-1) = fY;
            }
        }
        else{
            for (i=1;i<=runs;i++) {
                PT = PsiDraws.slice(i-1);
                BT = BetaDraws.slice(i-1);
                ST = SigmaDraws.slice(i-1);
                //
                fY.zeros();
                kY = Y0;
                //
                for (j=1;j<=periods;j++) {
                    fY.row(j-1) = kY*BT + PT - dX*(arma::kron(arma::eye(p,p),PT))*BT;
                    if (K > M){
                        kY(0,arma::span(M,K-1)) = kY(0,arma::span(0,K-M-1));
                    }
                    kY(0,arma::span(0,M-1)) = fY.row(j-1);
                }
                Forecasts.slice(i-1) = fY;
            }
        }
        //
        return Rcpp::List::create(Rcpp::Named("Forecasts") = Forecasts);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARS Forecast C++ exception (unknown reason)" );
    }
    return R_NilValue;
}