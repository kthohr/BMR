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

#include "bvarwc.h"
using namespace Rcpp;

SEXP WBVARReps(SEXP mSigma, SEXP mX, SEXP mZ, SEXP mY, SEXP maPr, SEXP mSPr, SEXP mvPr,
               SEXP mBVPr, SEXP mTp, SEXP mM, SEXP mK, SEXP mburnin, SEXP mkreps)
{
    try {
        arma::mat Sigma = as<arma::mat>(mSigma);
        arma::mat X = as<arma::mat>(mX);
        arma::mat Z = as<arma::mat>(mZ);
        arma::mat Y = as<arma::mat>(mY);
        arma::mat aPr = as<arma::mat>(maPr);
        arma::mat SPr = as<arma::mat>(mSPr);
        int vPr = as<int>(mvPr);
        arma::mat BVPr = as<arma::mat>(mBVPr);
        //
        int Tp = as<int>(mTp);
        int M = as<int>(mM);
        int K = as<int>(mK);
        int burnin = as<int>(mburnin);
        int kreps = as<int>(mkreps);
        //
        arma::cube BetaDraws(K, M, kreps);
        BetaDraws.zeros();
        arma::cube SigmaDraws(M, M, kreps);
        SigmaDraws.zeros();
        //
        arma::mat invBVPr = arma::inv(BVPr);
        arma::colvec YStacked = arma::vectorise(Y);
        int vPt = Tp + vPr;
        //
        int i,j;
        //
        arma::mat SigmaD = arma::kron(arma::inv_sympd(Sigma),arma::eye(Tp,Tp));
        arma::mat VPt = invBVPr + Z.t()*SigmaD*Z;
        //
        arma::mat VQ, VR;
        arma::qr(VQ,VR,VPt);
        VPt = arma::inv(VR)*VQ.t();
        //
        arma::mat aPt = VPt*(invBVPr*aPr + Z.t()*SigmaD*YStacked);
        arma::vec Krand = arma::randn(M*K);
        arma::vec alpha = aPt + arma::trans(arma::chol(VPt))*Krand;
        arma::mat Beta = arma::reshape(alpha,K,M);
        //
        arma::mat epsilon = Y - X*Beta;
        arma::mat SPt = SPr + epsilon.t()*epsilon;
        arma::mat Krand2 = arma::randn(SPt.n_rows,vPt);
        arma::mat A = arma::trans(arma::chol(inv_sympd(SPt)))*Krand2;
        Sigma = arma::inv_sympd(A*A.t());
        //
        for(i=1; i<=burnin; i++){
            SigmaD = arma::kron(inv_sympd(Sigma),arma::eye(Tp,Tp));
            VPt = invBVPr + Z.t()*SigmaD*Z;
            //
            VQ.zeros(), VR.zeros();
            arma::qr(VQ,VR,VPt);
            VPt = arma::inv(VR)*VQ.t();
            //
            aPt = VPt*(invBVPr*aPr + Z.t()*SigmaD*YStacked);
            Krand = arma::randn(M*K);
            alpha = aPt + arma::trans(arma::chol(VPt))*Krand;
            Beta = arma::reshape(alpha,K,M);
            //
            epsilon = Y - X*Beta;
            SPt = SPr + epsilon.t()*epsilon;
            Krand2 = arma::randn(SPt.n_rows,vPt);
            A = arma::trans(arma::chol(inv_sympd(SPt)))*Krand2;
            Sigma = arma::inv_sympd(A*A.t());
        }
        //
        for(j=1; j<=kreps; j++){
            SigmaD = arma::kron(inv_sympd(Sigma),arma::eye(Tp,Tp));
            VPt = invBVPr + Z.t()*SigmaD*Z;
            //
            VQ.zeros(), VR.zeros();
            arma::qr(VQ,VR,VPt);
            VPt = arma::inv(VR)*VQ.t();
            //
            aPt = VPt*(invBVPr*aPr + Z.t()*SigmaD*YStacked);
            Krand = arma::randn(M*K);
            alpha = aPt + arma::trans(arma::chol(VPt))*Krand;
            Beta = arma::reshape(alpha,K,M);
            //
            epsilon = Y - X*Beta;
            SPt = SPr + epsilon.t()*epsilon;
            Krand2 = arma::randn(SPt.n_rows,vPt);
            A = arma::trans(arma::chol(inv_sympd(SPt)))*Krand2;
            Sigma = arma::inv_sympd(A*A.t());
            //
            BetaDraws.slice(j-1) = Beta;
            SigmaDraws.slice(j-1) = Sigma;
        }
        //
        return Rcpp::List::create(Rcpp::Named("Beta") = BetaDraws,Rcpp::Named("Sigma") = SigmaDraws);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARW Gibbs C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP WBVARRepsB(SEXP mSigma, SEXP mX, SEXP mZ, SEXP mY, SEXP maPr, SEXP mSPr, SEXP mvPr, 
                SEXP mBVPr, SEXP mTp, SEXP mM, SEXP mK, SEXP mburnin)
{
    try {
        arma::mat Sigma = as<arma::mat>(mSigma);
        arma::mat X = as<arma::mat>(mX);
        arma::mat Z = as<arma::mat>(mZ);
        arma::mat Y = as<arma::mat>(mY);
        arma::mat aPr = as<arma::mat>(maPr);
        arma::mat SPr = as<arma::mat>(mSPr);
        int vPr = as<int>(mvPr);
        arma::mat BVPr = as<arma::mat>(mBVPr);
        //
        int Tp = as<int>(mTp);
        int M = as<int>(mM);
        int K = as<int>(mK);
        int burnin = as<int>(mburnin);
        //
        arma::mat invBVPr = arma::inv(BVPr);
        arma::colvec YStacked = arma::vectorise(Y);
        int vPt = Tp + vPr;
        //
        int i;
        //
        arma::mat SigmaD = arma::kron(inv_sympd(Sigma),arma::eye(Tp,Tp));
        arma::mat VPt = invBVPr + Z.t()*SigmaD*Z;
        //
        arma::mat VQ, VR;
        arma::qr(VQ,VR,VPt);
        VPt = arma::inv(VR)*VQ.t();
        //
        arma::mat aPt = VPt*(invBVPr*aPr + Z.t()*SigmaD*YStacked);
        arma::vec Krand = arma::randn(M*K);
        arma::vec alpha = aPt + arma::trans(arma::chol(VPt))*Krand;
        arma::mat Beta = arma::reshape(alpha,K,M);
        //
        arma::mat epsilon = Y - X*Beta;
        arma::mat SPt = SPr + epsilon.t()*epsilon;
        arma::mat Krand2 = arma::randn(SPt.n_rows,vPt);
        arma::mat A = arma::trans(arma::chol(inv_sympd(SPt)))*Krand2;
        Sigma = arma::inv_sympd(A*A.t());
        //
        for(i=1;i<=burnin;i++){
            SigmaD = arma::kron(inv_sympd(Sigma),arma::eye(Tp,Tp));
            VPt = invBVPr + Z.t()*SigmaD*Z;
            //
            VQ.zeros(), VR.zeros();
            arma::qr(VQ,VR,VPt);
            VPt = arma::inv(VR)*VQ.t();
            //
            aPt = VPt*(invBVPr*aPr + Z.t()*SigmaD*YStacked);
            Krand = arma::randn(M*K);
            alpha = aPt + arma::trans(arma::chol(VPt))*Krand;
            Beta = arma::reshape(alpha,K,M);
            //
            epsilon = Y - X*Beta;
            SPt = SPr + epsilon.t()*epsilon;
            Krand2 = arma::randn(SPt.n_rows,vPt);
            A = arma::trans(arma::chol(inv_sympd(SPt)))*Krand2;
            Sigma = arma::inv_sympd(A*A.t());
        }
        //
        return Rcpp::List::create(Rcpp::Named("Beta") = Beta,Rcpp::Named("Sigma") = Sigma);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARW Gibbs C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP WBVARRepsK(SEXP mSigma, SEXP mX, SEXP mZ, SEXP mY, SEXP maPr, SEXP mSPr, SEXP mvPr, 
                SEXP mBVPr, SEXP mTp, SEXP mM, SEXP mK, SEXP mkreps)
{
    try {
        arma::mat Sigma = as<arma::mat>(mSigma);
        arma::mat X = as<arma::mat>(mX);
        arma::mat Z = as<arma::mat>(mZ);
        arma::mat Y = as<arma::mat>(mY);
        arma::mat aPr = as<arma::mat>(maPr);
        arma::mat SPr = as<arma::mat>(mSPr);
        int vPr = as<int>(mvPr);
        arma::mat BVPr = as<arma::mat>(mBVPr);
        //
        int Tp = as<int>(mTp);
        int M = as<int>(mM);
        int K = as<int>(mK);
        int kreps = as<int>(mkreps);
        //
        arma::cube BetaDraws(K, M, kreps);
        BetaDraws.zeros();
        arma::cube SigmaDraws(M, M, kreps);
        SigmaDraws.zeros();
        //
        arma::mat invBVPr = arma::inv(BVPr);
        arma::colvec YStacked = arma::vectorise(Y);
        int vPt = Tp + vPr;
        //
        int j;
        //
        arma::mat SigmaD = arma::kron(inv_sympd(Sigma),arma::eye(Tp,Tp));
        arma::mat VPt = invBVPr + Z.t()*SigmaD*Z;
        //
        arma::mat VQ, VR;
        arma::qr(VQ,VR,VPt);
        VPt = arma::inv(VR)*VQ.t();
        //
        arma::mat aPt = VPt*(invBVPr*aPr + Z.t()*SigmaD*YStacked);
        arma::vec Krand = arma::randn(M*K);
        arma::vec alpha = aPt + arma::trans(arma::chol(VPt))*Krand;
        arma::mat Beta = arma::reshape(alpha,K,M);
        //
        arma::mat epsilon = Y - X*Beta;
        arma::mat SPt = SPr + epsilon.t()*epsilon;
        arma::mat Krand2 = arma::randn(SPt.n_rows,vPt);
        arma::mat A = arma::trans(arma::chol(inv_sympd(SPt)))*Krand2;
        Sigma = arma::inv_sympd(A*A.t());
        //
        for (j=1;j<=kreps;j++) {
            SigmaD = arma::kron(inv_sympd(Sigma),arma::eye(Tp,Tp));
            VPt = invBVPr + Z.t()*SigmaD*Z;
            //
            VQ.zeros(), VR.zeros();
            arma::qr(VQ,VR,VPt);
            VPt = arma::inv(VR)*VQ.t();
            //
            aPt = VPt*(invBVPr*aPr + Z.t()*SigmaD*YStacked);
            Krand = arma::randn(M*K);
            alpha = aPt + arma::trans(arma::chol(VPt))*Krand;
            Beta = arma::reshape(alpha,K,M);
            //
            epsilon = Y - X*Beta;
            SPt = SPr + epsilon.t()*epsilon;
            Krand2 = arma::randn(SPt.n_rows,vPt);
            A = arma::trans(arma::chol(inv_sympd(SPt)))*Krand2;
            Sigma = arma::inv_sympd(A*A.t());
            //
            BetaDraws.slice(j-1) = Beta;
            SigmaDraws.slice(j-1) = Sigma;
        }
        //
        return Rcpp::List::create(Rcpp::Named("Beta") = BetaDraws,Rcpp::Named("Sigma") = SigmaDraws);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARW Gibbs C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP WBVARIRFs(SEXP mM, SEXP mK, SEXP mkcons, SEXP mkreps, SEXP mirfperiods,
               SEXP mBDs, SEXP mSDs)
{
    try {
        int M = as<int>(mM);
        int K = as<int>(mK);
        int kcons = as<int>(mkcons);
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
        arma::mat BetaT = arma::zeros<arma::mat>(K-kcons,M);
        arma::mat shock = arma::zeros<arma::mat>(M,M);
        arma::mat ShockO = arma::zeros<arma::mat>(K-M-kcons,M);
        arma::mat ShockB = arma::zeros<arma::mat>(K-kcons,M);
        arma::mat Pow = arma::zeros<arma::mat>(M,M);
        //
        for(j=1; j<=kreps; j++){
            BetaT = BetaDraws(arma::span(kcons,K-1),arma::span(),arma::span(j-1,j-1));
            BetaT = BetaT.t();
            //
            shock = SigmaDraws.slice(j-1);
            shock = arma::trans(arma::chol(shock));
            //
            ShockB.zeros();
            ShockB(arma::span(0,M-1),arma::span()) = shock;
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
        ::Rf_error( "BMR: BVARW IRF C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP bvarwforecast(SEXP mY0, SEXP mM, SEXP mK, SEXP mkcons, SEXP mruns,
                   SEXP mperiods, SEXP minclshocks, SEXP mBDs, SEXP mSDs)
{
    try {
        int M = as<int>(mM);
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
        arma::mat BT = BetaDraws.slice(0);
        arma::mat ST = SigmaDraws.slice(0);
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
            for (i=1;i<=runs;i++) {
                BT = BetaDraws.slice(i-1);
                ST = SigmaDraws.slice(i-1);
                //
                Shock = arma::trans(arma::chol(ST));
                //
                fY.zeros();
                kY = Y0;
                //
                for (j=1;j<=periods;j++) {
                    Krand = arma::randn(M);
                    jshock = arma::trans(Shock*Krand);
                    //
                    fY.row(j-1) = kY*BT + jshock;
                    if (K > M+kcons){
                        kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1));
                    }
                    kY(0,arma::span(kcons,M-1+kcons)) = fY.row(j-1);
                }
                Forecasts.slice(i-1) = fY;
            }
        }
        else{
            for (i=1;i<=runs;i++) {
                BT = BetaDraws.slice(i-1);
                ST = SigmaDraws.slice(i-1);
                //
                fY.zeros();
                kY = Y0;
                //
                for (j=1;j<=periods;j++) {
                    fY.row(j-1) = kY*BT;
                    if (K > M+kcons){
                        kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1));
                    }
                    kY(0,arma::span(kcons,M-1+kcons)) = fY.row(j-1);
                }
                Forecasts.slice(i-1) = fY;
            }
        }
        //
        return Rcpp::List::create(Rcpp::Named("Forecasts") = Forecasts);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARW IRF C++ exception (unknown reason)" );
    }
    return R_NilValue;
}