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

#include "dsgevarc.h"
using namespace Rcpp;

SEXP DSGEVARPriorC(SEXP mdsgedata, SEXP mObserveMat, SEXP mObsCons, SEXP mF, SEXP mG,
                   SEXP mshocks, SEXP mR, SEXP mp, SEXP mmax_iter)
{
    try {
        arma::mat dsgedata = as<arma::mat>(mdsgedata);
        //
        arma::mat ObserveMat = as<arma::mat>(mObserveMat);
        arma::mat ObsCons = as<arma::mat>(mObsCons);
        //
        arma::mat F = as<arma::mat>(mF);
        arma::mat G = as<arma::mat>(mG);
        arma::mat Q = as<arma::mat>(mshocks);
        arma::mat R = as<arma::mat>(mR);
        //
        arma::mat OC2 = ObsCons*ObsCons.t();
        //
        int p = as<int>(mp);
        int max_iter = as<int>(mmax_iter);
        //
        int i,j;
        //
        arma::mat GQG = G*Q*G.t();
        //
        arma::mat SigmaSS = arma::eye(F.n_cols,F.n_cols);
        //
        arma::mat SigmaSSOld = GQG;
        arma::mat SigmaSSNew = GQG;
        arma::mat IMat = F;
        //
        for(j=1;j<=max_iter;j++){
            SigmaSSNew = SigmaSSOld + IMat*SigmaSSOld*IMat.t();
            IMat *= IMat;
            //
            SigmaSSOld = SigmaSSNew;
        }
        SigmaSS = SigmaSSNew;
        //
        arma::cube SigmaX(dsgedata.n_cols, dsgedata.n_cols, p+1);
        SigmaX.zeros();
        //
        SigmaX.slice(0) = OC2 + ObserveMat.t()*SigmaSS*ObserveMat + R;
        //
        arma::mat Fp = arma::eye(F.n_cols,F.n_cols);
        for(i=1;i<=p;i++){
            Fp *= F;
            SigmaX.slice(i) = OC2 + ObserveMat.t()*Fp*SigmaSS*ObserveMat;
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

SEXP DSGEVARLikelihood(SEXP mlogGPR, SEXP mXX, SEXP mGammaYY, SEXP mGammaXY, SEXP mGammaXX,
                       SEXP mGammaBarYY, SEXP mGammaBarXY, SEXP mGammaBarXX,
                       SEXP mlambda, SEXP mTp, SEXP mM, SEXP mp, SEXP mkcons)
{
    try {
        //
        double logGPR = as<double>(mlogGPR);
        double lambda = as<double>(mlambda);
        int Tp = as<int>(mTp);
        int M = as<int>(mM);
        int p = as<int>(mp);
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
        double lambdaT = Tp*lambda;
        double invlambda = 1/lambda;
        double M2 = M/2;
        //
        arma::mat invGammaXX = arma::inv_sympd(GammaXX);
        arma::mat invGammaBarXX = arma::inv_sympd(GammaBarXX);
        arma::mat SigmaEpsilon = GammaYY - GammaXY*invGammaXX*GammaXY.t();
        arma::mat SigmaHatEpsilon = GammaBarYY - GammaBarXY.t()*invGammaBarXX*GammaBarXY;
        //
        arma::mat invlambdaXX = invlambda*XX;
        arma::mat invlambdaSigma = (1 + invlambda)*SigmaHatEpsilon;
        //
        double Term1 = -M2*log(arma::det(GammaXX + invlambdaXX)) + M2*log(arma::det(GammaXX));
        double Term2 = -((Tp + lambdaT - M*p - kcons)/2)*log(arma::det( invlambdaSigma ));
        double Term3 = ((lambdaT - M*p - kcons)/2)*log(arma::det(SigmaEpsilon));
        //
        double logLikelihood = Term1 + Term2 + Term3 + logGPR - (((M*Tp)/2)*log(lambdaT*arma::datum::pi));
        // (note that this is not the negative of the loglikelihood)
        return Rcpp::List::create(Rcpp::Named("logLikelihood") = logLikelihood);
    } catch( std::exception &ex ) {
        //forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: DSGEVAR Likelihood C++ exception (unknown reason)" );
    }
    //return R_NilValue;
    double LLKerr = -1000000;
    return Rcpp::List::create(Rcpp::Named("logLikelihood") = LLKerr);
}

SEXP DSGEVARLikelihoodInf(SEXP mYY, SEXP mXY, SEXP mXX,
                          SEXP mGammaYY, SEXP mGammaXY, SEXP mGammaXX,
                          SEXP mTp, SEXP mM)
{
    try {
        //
        int Tp = as<int>(mTp);
        int M = as<int>(mM);
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
        arma::mat SigmaEpsilon = GammaYY - GammaXY*invGammaXX*GammaXY.t();
        arma::mat Beta = invGammaXX*GammaXY.t();
        arma::mat SigmaBarEpsilon = YY + Beta.t()*XX*Beta - Beta.t()*XY - XY.t()*Beta;
        //
        double logLikelihood = - (((M*Tp)/2)*log(2*arma::datum::pi)) - ((Tp/2)*log(arma::det(SigmaEpsilon))) - ((Tp/2)*arma::trace(arma::inv_sympd(SigmaEpsilon)*SigmaBarEpsilon));
        // (note that this is not the negative of the loglikelihood)
        return Rcpp::List::create(Rcpp::Named("logLikelihood") = logLikelihood);
    } catch( std::exception &ex ) {
        //forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: DSGEVAR Likelihood C++ exception (unknown reason)" );
    }
    //return R_NilValue;
    double LLKerr = -1000000;
    return Rcpp::List::create(Rcpp::Named("logLikelihood") = LLKerr);
}

SEXP DSGEVARReps(SEXP mGammaBarYY, SEXP mGammaBarXY, SEXP mGammaBarXX, SEXP mGXX,
                 SEXP mGHXX, SEXP mlambda, SEXP mkreps, SEXP mTp, SEXP mM, SEXP mp,
                 SEXP mkcons)
{
    try {
        //
        double lambda = as<double>(mlambda);
        int kreps = as<int>(mkreps);
        int Tp = as<int>(mTp);
        int M = as<int>(mM);
        int p = as<int>(mp);
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
        arma::mat GBarYY = GammaBarYY.slice(0);
        arma::mat GBarXY = GammaBarXY.slice(0);
        arma::mat GBarXX = GammaBarXX.slice(0);
        //
        arma::mat invGBarXX = arma::inv_sympd(GBarXX);
        arma::mat BetaHat = invGBarXX*GBarXY;
        arma::mat SigmaEpsilon = GBarYY - GBarXY.t()*invGBarXX*GBarXY;
        //
        double S2 = (1+lambda)*Tp;
        int S3 = (1+lambda)*Tp - M*p - kcons;
        //
        arma::mat SigmaDraw = S2*SigmaEpsilon;
        arma::mat KSigmaRand = arma::randn(SigmaDraw.n_rows,S3);
        arma::mat A = arma::trans(arma::chol(arma::inv_sympd(SigmaDraw)))*KSigmaRand;
        SigmaDraw = arma::inv_sympd(A*A.t());
        //
        double B1 = 1/lambda;
        double B2 = 1/(Tp*lambda);
        arma::mat GammaXX = GXX.slice(0);
        //
        arma::mat BetaCov = arma::inv_sympd(GammaXX + B1*GHXX);
        arma::vec KBetaRand = arma::randn(M*(M*p+kcons));
        arma::colvec VecBetaHat = arma::vectorise(BetaHat);
        arma::mat VecBetaDraw = VecBetaHat + arma::trans(arma::chol(B2*arma::kron(SigmaEpsilon,BetaCov)))*KBetaRand;
        arma::mat BetaDraw = arma::reshape(VecBetaDraw, M*p+kcons, M);
        //
        BetaDraws.slice(0) = BetaDraw;
        SigmaDraws.slice(0) = SigmaDraw;
        //
        for(j=2;j<=kreps;j++){
            GBarYY = GammaBarYY.slice(j-1);
            GBarXY = GammaBarXY.slice(j-1);
            GBarXX = GammaBarXX.slice(j-1);
            //
            invGBarXX = arma::inv_sympd(GBarXX);
            BetaHat = invGBarXX*GBarXY;
            SigmaEpsilon = GBarYY - GBarXY.t()*invGBarXX*GBarXY;
            //
            SigmaDraw = S2*SigmaEpsilon;
            KSigmaRand = arma::randn(SigmaDraw.n_rows,S3);
            A = arma::trans(arma::chol(arma::inv(SigmaDraw)))*KSigmaRand;
            SigmaDraw = arma::inv_sympd(A*A.t());
            //
            GammaXX = GXX.slice(j-1);
            //
            BetaCov = arma::inv_sympd(GammaXX + B1*GHXX);
            KBetaRand = arma::randn(M*(M*p+kcons));
            VecBetaHat = arma::vectorise(BetaHat);
            VecBetaDraw = VecBetaHat + arma::trans(arma::chol(B2*arma::kron(SigmaEpsilon,BetaCov)))*KBetaRand;
            BetaDraw = arma::reshape(VecBetaDraw, M*p+kcons, M);
            //
            BetaDraws.slice(j-1) = BetaDraw;
            SigmaDraws.slice(j-1) = SigmaDraw;
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

SEXP DSGEVARRepsInf(SEXP mGammaBarYY, SEXP mGammaBarXY, SEXP mGammaBarXX,
                    SEXP mkreps, SEXP mM, SEXP mp, SEXP mkcons)
{
    try {
        //
        int kreps = as<int>(mkreps);
        int M = as<int>(mM);
        int p = as<int>(mp);
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
        arma::mat GBarYY = GammaBarYY.slice(0);
        arma::mat GBarXY = GammaBarXY.slice(0);
        arma::mat GBarXX = GammaBarXX.slice(0);
        //
        arma::mat invGBarXX = arma::inv_sympd(GBarXX);
        arma::mat BetaHat = invGBarXX*GBarXY;
        arma::mat SigmaEpsilon = GBarYY - GBarXY.t()*invGBarXX*GBarXY;
        //
        BetaDraws.slice(0) = BetaHat;
        SigmaDraws.slice(0) = SigmaEpsilon;
        //
        for (j=2;j<=kreps;j++) {
            GBarYY = GammaBarYY.slice(j-1);
            GBarXY = GammaBarXY.slice(j-1);
            GBarXX = GammaBarXX.slice(j-1);
            //
            invGBarXX = arma::inv_sympd(GBarXX);
            BetaHat = invGBarXX*GBarXY;
            SigmaEpsilon = GBarYY - GBarXY.t()*invGBarXX*GBarXY;
            //
            BetaDraws.slice(j-1) = BetaHat;
            SigmaDraws.slice(j-1) = SigmaEpsilon;
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
        arma::mat BetaT = arma::zeros<arma::mat>(K-kcons,M);
        arma::mat shock = arma::zeros<arma::mat>(M,M);
        arma::mat ShockO = arma::zeros<arma::mat>(K-M-kcons,M);
        arma::mat ShockB = arma::zeros<arma::mat>(K-kcons,M);
        arma::mat Pow = arma::zeros<arma::mat>(M,M);
        //
        for(j=1;j<=kreps;j++){
            BetaT = BetaDraws(arma::span(kcons,K-1),arma::span(),arma::span(j-1,j-1));
            BetaT = BetaT.t();
            shock = SigmaDraws.slice(j-1);
            ShockB.zeros();
            ShockB.rows(0,M-1) = shock;
            ImpStore.slice((j-1)*irfperiods) = shock;
            for(i=2;i<=irfperiods;i++){
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
        ::Rf_error( "BMR: DSGEVAR IRF C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP dsgevarforecast(SEXP mY0, SEXP mM, SEXP mK, SEXP mkcons,
                     SEXP mruns, SEXP mperiods, SEXP minclshocks, SEXP mBDs, SEXP mSDs)
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
            for(i=1;i<=runs;i++){
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
            for (i=1;i<=runs;i++) {
                BT = BetaDraws.slice(i-1);
                ST = SigmaDraws.slice(i-1);
                //
                fY.zeros();
                kY = Y0;
                //
                for(j=1; j<=periods; j++){
                    fY.row(j-1) = kY*BT;
                    if(K > M+kcons){
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
        ::Rf_error( "BMR: DSGEVAR Forecast C++ exception (unknown reason)" );
    }
    return R_NilValue;
}