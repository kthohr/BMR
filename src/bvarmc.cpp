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

#include "bvarmc.h"
using namespace Rcpp;

SEXP MBVARReps(SEXP mSigma, SEXP mZ, SEXP mY, SEXP maPr, SEXP mBVPr, SEXP mM, SEXP mK,
               SEXP mburnin, SEXP mkreps)
{
    try {
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
        arma::colvec YStacked = arma::vectorise(Y);
        arma::mat invBVPr = arma::inv(BVPr);
        //
        arma::mat BVPt = arma::inv_sympd(invBVPr + Z.t()*Sigma*Z);
        arma::mat aPt = BVPt*(invBVPr*aPr + Z.t()*Sigma*YStacked);
        //
        arma::vec Krand = arma::randn(M*K);
        arma::mat tCholBVPt = arma::trans(arma::chol(BVPt));
        arma::vec alpha = aPt + tCholBVPt*Krand;
        //
        for(i=1; i<=burnin; i++){
            Krand = arma::randn(M*K);
            alpha = aPt + tCholBVPt*Krand;
        }
        //
        for(j=1; j<=kreps; j++){
            Krand = arma::randn(M*K);
            alpha = aPt + tCholBVPt*Krand;
            //
            BetaDraws.slice(j-1) = arma::reshape(alpha,K,M);
        }
        //
        return Rcpp::List::create(Rcpp::Named("Beta") = BetaDraws);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: BVARM Gibbs C++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP MBVARIRFs( SEXP mshock, SEXP mM, SEXP mK, SEXP mkcons, SEXP mkreps, 
                SEXP mirfperiods, SEXP mBDs )
{
    try {
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
        shock = shock.t();
        //
        arma::mat BetaT = arma::zeros<arma::mat>(K-kcons,M);
        arma::mat ShockO = arma::zeros<arma::mat>(K-M-kcons,M);
        arma::mat ShockB = arma::zeros<arma::mat>(K-kcons,M);
        arma::mat Pow = arma::zeros<arma::mat>(M,M);
        //
        ShockB(arma::span(0,M-1),arma::span()) = shock;
        //
        for(j=1; j<=kreps; j++){
            BetaT = BetaDraws(arma::span(kcons,K-1),arma::span(),arma::span(j-1,j-1));
            BetaT = BetaT.t();
            ShockB.zeros();
            ShockB.rows(0,M-1) = shock;
            ImpStore.slice((j-1)*irfperiods) = shock;
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
        ::Rf_error( "BMR: BVARM IRF C++ exception (unknown reason)" );
    }
	return R_NilValue;
}

SEXP bvarmforecast(SEXP mY0, SEXP mM, SEXP mK, SEXP mkcons, SEXP mruns, SEXP mperiods,
                   SEXP minclshocks, SEXP mBDs, SEXP mShock)
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
        arma::cube BetaDraws(BDArray.begin(), K, M, runs, false);
        //
        arma::mat Shock = as<arma::mat>(mShock);
        //
        arma::cube Forecasts(periods, M, runs);
        Forecasts.zeros();
        //
        arma::mat BT = BetaDraws.slice(0);
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
                BT = BetaDraws.slice(i-1);
                fY.zeros();
                kY = Y0;
                for (j=1;j<=periods;j++) {
                    Krand = arma::randn(M);
                    jshock = arma::trans(Shock*Krand);
                    //
                    fY.row(j-1) = kY*BT + jshock;
                    if(K > M+kcons){
                        kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1));
                    }
                    kY(0,arma::span(kcons,M-1+kcons)) = fY.row(j-1);
                }
                Forecasts.slice(i-1) = fY;
            }
        }
        else{
            for(i=1; i<=runs; i++){
                BT = BetaDraws.slice(i-1);
                fY.zeros();
                kY = Y0;
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
        ::Rf_error( "BMR: BVARM Forecast C++ exception (unknown reason)" );
    }
    return R_NilValue;
}