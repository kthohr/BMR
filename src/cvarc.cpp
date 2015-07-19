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

#include "cvarc.h"
using namespace Rcpp;

SEXP CVARReps(SEXP mBeta, SEXP mSigma, SEXP mX, SEXP mY, SEXP mTp, SEXP mM,
              SEXP mK, SEXP mkcons, SEXP mboot, SEXP mRDs)
{
    try {
        arma::mat Y = as<arma::mat>(mY);
        arma::mat X = as<arma::mat>(mX);
        arma::mat Beta = as<arma::mat>(mBeta);
        arma::mat Sigma = as<arma::mat>(mSigma);
        //
        double Tp = as<int>(mTp);
        int M = as<int>(mM);
        double K = as<int>(mK);
        int kcons = as<int>(mkcons);
        int boot = as<int>(mboot);
        //
        NumericVector RDArray(mRDs);
        arma::cube ResidDraws(RDArray.begin(), Tp, M, boot, false);
        //
        arma::cube BetaDraws(K, M, boot);
        BetaDraws.zeros();
        arma::cube SigmaDraws(M, M, boot);
        SigmaDraws.zeros();
        //
        double SigmaS2 = 1/(Tp-K);
        //
        arma::mat BetaS = Beta;
        arma::mat SigmaS = Sigma;
        arma::mat Epsilon = Y;
        Epsilon.zeros();
        //
        int i,j;
        //
        arma::mat YT = arma::zeros<arma::mat>(Tp,M);
        arma::mat XT = arma::zeros<arma::mat>(Tp,K);
        //
        arma::mat XTM = arma::zeros<arma::mat>(1,K);
        XTM(0,0) = 1;
        //
        arma::mat RT = YT;
        //
        for(i=1; i<=boot; i++){
            RT = ResidDraws.slice(i-1);
            XTM = X.row(0);
            //
            for(j=1; j<=Tp; j++){
                XT.row(j-1) = XTM;
                YT.row(j-1) = XTM*Beta + RT.row(j-1);
                //
                if(K > M+kcons){
                    XTM(0,arma::span(M+kcons,K-1)) = XTM(0,arma::span(kcons,K-M-1));
                }
                XTM(0,arma::span(kcons,M-1+kcons)) = YT.row(j-1);
            }
            //
            //BetaS = arma::solve(XT.t()*XT,XT.t()*YT);
            BetaS = arma::solve(XT,YT);
            //
            Epsilon = YT - XT*BetaS;
            SigmaS = Epsilon.t()*Epsilon;
            SigmaS = SigmaS2*SigmaS;
            //
            BetaDraws.slice(i-1) = BetaS;
            SigmaDraws.slice(i-1) = SigmaS;
        }
        //
        return Rcpp::List::create(Rcpp::Named("Beta") = BetaDraws,Rcpp::Named("Sigma") = SigmaDraws);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: CVAR Bootstrap C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP CVARIRFs(SEXP mM, SEXP mK, SEXP mkcons, SEXP mkreps, SEXP mirfperiods,
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
        arma::cube BetaDraws(BDArray.begin(), K, M, kreps, false);
        arma::cube SigmaDraws(SDArray.begin(), M, M, kreps, false);
        //
        arma::cube ImpStore(M, M, irfperiods*kreps);
        ImpStore.zeros();
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
        ::Rf_error( "BMR: CVAR IRF C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP cvarforecast(SEXP mY0, SEXP mM, SEXP mp, SEXP mK, SEXP mkcons, SEXP mperiods,
                  SEXP mci, SEXP mBeta, SEXP mSigma)
{
    try {
        int M = as<int>(mM);
        int p = as<int>(mp);
        int K = as<int>(mK);
        int kcons = as<int>(mkcons);
        int periods = as<int>(mperiods);
        double ci = as<double>(mci);
        //
        arma::mat Y0 = as<arma::mat>(mY0);
        //
        arma::mat Beta = as<arma::mat>(mBeta);
        arma::mat BetaNC = Beta.rows(kcons,K-1);
        arma::mat Sigma = as<arma::mat>(mSigma);
        //
        arma::cube Forecasts(periods, M, 3);
        Forecasts.zeros();
        //
        arma::mat fY = arma::zeros<arma::mat>(periods,M);
        arma::mat fYL = arma::zeros<arma::mat>(periods,M);
        arma::mat fYU = arma::zeros<arma::mat>(periods,M);
        arma::mat kY = Y0;
        //
        arma::vec SEY = arma::zeros<arma::mat>(M,1);
        //
        int j;
        //
        arma::mat Shocks = arma::zeros<arma::mat>(p*M,p*M);
        arma::mat SEForecast = arma::zeros<arma::mat>(M,M);
        //
        for(j=1; j<=periods; j++){
            SEForecast = Sigma + BetaNC.t()*Shocks*BetaNC;
            if(K > M+kcons){
                Shocks(arma::span(M,(M*p)-1),arma::span(M,(M*p)-1)) = Shocks(arma::span(0,(M*(p-1))-1),arma::span(0,(M*(p-1))-1));
            }
            Shocks(arma::span(0,M-1),arma::span(0,M-1)) = Sigma;
            //
            SEY = SEForecast.diag();
            //
            fY.row(j-1) = kY*Beta;
            fYL.row(j-1) = kY*Beta - ci*arma::trans(sqrt(SEY));
            fYU.row(j-1) = kY*Beta + ci*arma::trans(sqrt(SEY));
            //
            if(K > M+kcons){
                kY(0,arma::span(M+kcons,K-1)) = kY(0,arma::span(kcons,K-M-1));
            }
            //
            kY(0,arma::span(kcons,M-1+kcons)) = fY.row(j-1);
        }
        // Point Forecast, then lower CI, then upper CI
        Forecasts.slice(0) = fY;
        Forecasts.slice(1) = fYL;
        Forecasts.slice(2) = fYU;
        //
        return Rcpp::List::create(Rcpp::Named("Forecasts") = Forecasts,Rcpp::Named("Shocks") = Shocks);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: CVAR Forecast C++ exception (unknown reason)" );
    }
    return R_NilValue;
}