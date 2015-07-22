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

#include "dsgec.h"
using namespace Rcpp;

#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
# define PREDEF_PLATFORM_UNIX
#endif

SEXP UhligCpp(SEXP mA, SEXP mB, SEXP mC, SEXP mD,
              SEXP mF, SEXP mG, SEXP mH,
              SEXP mJ, SEXP mK, SEXP mL, SEXP mM, SEXP mN,
              SEXP mwhichEig)
{
    try {
        //Suppress warnings:
        std::ostream nullstream(0);
        arma::set_stream_err2(nullstream);
        //
        arma::mat A = as<arma::mat>(mA);
        arma::mat B = as<arma::mat>(mB);
        arma::mat C = as<arma::mat>(mC);
        arma::mat D = as<arma::mat>(mD);
        arma::mat F = as<arma::mat>(mF);
        arma::mat G = as<arma::mat>(mG);
        arma::mat H = as<arma::mat>(mH);
        arma::mat J = as<arma::mat>(mJ);
        arma::mat K = as<arma::mat>(mK);
        arma::mat L = as<arma::mat>(mL);
        arma::mat M = as<arma::mat>(mM);
        arma::mat N = as<arma::mat>(mN);
        //
        arma::vec which_eig = as<arma::vec>(mwhichEig);
        //
        arma::mat P, Q, R, S;
        arma::cx_vec eigenvals; arma::cx_mat eigenvecs;
        //
        int l = C.n_rows;
        int n = C.n_cols;
        int k = std::min(N.n_rows,N.n_cols);
        //
        int m = 0;
        if( l == 0 ){
            m = F.n_cols;
        }else{ // l > 0
            m = A.n_cols;
        }
        //
        int pick_eig = 0;
        if(which_eig.n_elem > 0){
            pick_eig = 1;
        }
        //
        double bignum = 1e+07;
        //
        arma::mat Xi(2*m,2*m);
        arma::mat Delta(2*m,2*m);
        //
        arma::mat Psi = F;
        arma::mat Gamma = -G;
        arma::mat Theta = -H;
        //
        arma::mat Cpinv = C;
        if(l > 0){
            Cpinv = arma::pinv(C);
            Psi = F - J*Cpinv*A;
            Gamma = J*Cpinv*B - G + K*Cpinv*A;
            Theta = K*Cpinv*B - H;
        }
        //
        Xi(arma::span(0,m-1),arma::span(0,m-1)) = Gamma;
        Xi(arma::span(0,m-1),arma::span(m,2*m-1)) = Theta;
        Xi(arma::span(m,2*m-1),arma::span(0,m-1)) = arma::eye(m,m);
        Xi(arma::span(m,2*m-1),arma::span(m,2*m-1)).zeros();
        //
        Delta(arma::span(0,m-1),arma::span(0,m-1)) = Psi;
        Delta(arma::span(0,m-1),arma::span(m,2*m-1)).zeros();
        Delta(arma::span(m,2*m-1),arma::span(0,m-1)).zeros();
        Delta(arma::span(m,2*m-1),arma::span(m,2*m-1)) = arma::eye(m,m);
        /*
         * Next perform a generalized eigen decomp:
         */
        arma::cx_vec EigValue;
        arma::cx_mat EigVec;
        // generalized eigen decomp:
        arma::eig_pair(EigValue, EigVec, Xi, Delta);
        //
        arma::vec EigValueAbs = abs(EigValue);
        arma::vec EigValueReal = real(EigValue);
        /*
         * Deal with infinite values (otherwise sort will return an error):
         */
        arma::uvec infindices = find_nonfinite(EigValueAbs);
        arma::uvec neginfind = arma::find(EigValueReal < -bignum);
        //
        int infin2 = infindices.n_elem;
        int infin3 = neginfind.n_elem;
        //
        if(infin2 > 0){
            arma::cx_vec BigNum(infin2);
            BigNum.ones();
            //
            EigValueAbs.elem(infindices) = arma::ones(infin2)*bignum;
            EigValue.elem(infindices) = BigNum*bignum;
        }
        if(infin3 > 0){
            arma::cx_vec BigNum2(infin3);
            BigNum2.ones();
            //
            EigValue.elem(neginfind) = - BigNum2*bignum;
        }
        /*
         * Now sort the eigenvalues and eigenvectors...
         */
        arma::uvec indices = arma::sort_index(EigValueAbs);
        // ... from smallest to largest in absolute value
        arma::vec EigValueAbsSorted = EigValueAbs.elem(indices);
        arma::cx_vec EigValueSorted = EigValue.elem(indices);
        arma::cx_mat EigVecSorted = EigVec.cols(indices);
        /*
         * If the user prefers to 'choose' which eigenvalues to use...
         */
        if(pick_eig > 0){
            /*
             * by 'egvecind-1' this implies that the elements of 'which_eig' begin at 1 and not 0
             */
            indices = arma::conv_to<arma::uvec>::from(which_eig - 1);
            //
            EigValueSorted = EigValueSorted.elem(indices);
            EigVecSorted = EigVecSorted.cols(indices);
        }
        //
        eigenvals = EigValueSorted;
        eigenvecs = EigVecSorted;
        //
        arma::cx_vec LambdaVec = EigValueSorted.rows(0,m-1);
        arma::cx_mat Lambda = arma::diagmat(LambdaVec);
        //
        arma::cx_mat Omega = EigVecSorted(arma::span(m,2*m-1),arma::span(0,m-1));
        //
        P = arma::real( Omega*Lambda*arma::pinv(Omega) );
        /*
         * Now calculate Q, R, and S
         */
        Cpinv = C;
        //
        arma::mat LNM = L*N + M;
        arma::colvec LNMStacked = arma::vectorise(LNM);
        //
        arma::mat V = arma::kron(arma::trans(N),F) + arma::kron(arma::eye(k,k),(F*P+G));
        //
        arma::vec QS;
        R.set_size(0,P.n_cols); R.zeros();
        if(l == 0){
            QS = -arma::inv(V)*LNMStacked;
        }else{
            /*
             * First R, ...
             */
            Cpinv = pinv(C);
            R = - Cpinv*(A*P + B);
            /*
             * then Q
             */
            arma::mat V2 = arma::zeros<arma::mat>((k*A.n_rows) + V.n_rows, (k*A.n_rows) + V.n_rows);
            //
            V2(arma::span(0,(k*A.n_rows)-1),arma::span(0,(k*A.n_cols)-1)) = arma::kron(arma::eye(k,k),A);
            V2(arma::span(0,(k*A.n_rows)-1),arma::span(k*A.n_cols,(k*A.n_cols)+(k*C.n_cols)-1)) = arma::kron(arma::eye(k,k),C);
            V2(arma::span(k*A.n_rows,(k*A.n_rows)+V.n_rows-1),arma::span(0,(k*F.n_cols)-1)) = arma::kron(trans(N),F) + arma::kron(arma::eye(k,k),(F*P + J*R + G));
            V2(arma::span(k*A.n_rows,(k*A.n_rows)+V.n_rows-1),arma::span(k*F.n_cols,(k*F.n_cols)+(k*J.n_cols)-1)) = arma::kron(trans(N),J) + arma::kron(arma::eye(k,k),K);
            //
            arma::colvec DStacked = arma::vectorise(D);
            //
            arma::mat DLNM(DStacked.n_rows + LNMStacked.n_rows,1);
            DLNM.rows(0,DStacked.n_rows-1) = DStacked;
            DLNM.rows(DStacked.n_rows,DStacked.n_rows + LNMStacked.n_rows - 1) = LNMStacked;
            //
            QS = - arma::inv(V2)*DLNM;
        }
        arma::vec QVec = QS.rows(0,(m*k)-1);
        Q = arma::reshape(QVec, m, k);
        /*
         * Finally, S...
         */
        S.set_size(0,Q.n_cols); S.zeros();
        if(l > 0){
            arma::vec Sc = QS.rows(m*k,(m+n)*k-1);
            S = arma::reshape(Sc, n, k);
        }
        //
        return Rcpp::List::create(Rcpp::Named("P") = P, Rcpp::Named("Q") = Q, Rcpp::Named("R") = R,Rcpp::Named("S") = S, Rcpp::Named("EigValueSorted") = eigenvals, Rcpp::Named("EigVecSorted") = eigenvecs);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: SDSGE C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

int qzswitch(int i, arma::cx_mat &A, arma::cx_mat &B, arma::cx_mat &Q, arma::cx_mat &Z){
    //
    double realsmall = 1e-08;
    //
    std::complex<double> a = A(i,i);         std::complex<double> d = B(i,i);
    std::complex<double> b = A(i,i+1);     std::complex<double> e = B(i,i+1);
    std::complex<double> c = A(i+1,i+1); std::complex<double> f = B(i+1,i+1);
    //
    arma::cx_mat wz(2,1); arma::cx_mat wz1(2,1); arma::cx_mat wz2(2,2); arma::cx_mat wzt;
    arma::cx_mat xy(1,2); arma::cx_mat xy1; arma::cx_mat xy2(2,2); arma::cx_mat xyt;
    //
    std::complex<double> n,m;
    //
    //
    //
    if( std::abs(c) < realsmall && std::abs(f) < realsmall ){
        if( std::abs(a) < realsmall ){
            // l.r. coincident 0's with u.l. of A=0; do nothing
            return 0;
        }
        else{
            // l.r. coincident zeros; put 0 in u.l. of a
            wz.set_size(2,1);
            //
            wz(0,0) = b;
            wz(1,0) = -a;
            //
            wz = wz / arma::as_scalar(arma::sqrt(arma::trans(wz)*wz));
            wz1 = wz;
            wzt = arma::trans(wz);
            //
            wz2(arma::span(),0) = wz1;
            //
            wz2(0,1) = wzt(0,1);
            wz2(1,1) = -wzt(0,0);
            //
            xy2 = arma::eye<arma::cx_mat>(2,2);
        }
    }
    else if( std::abs(a) < realsmall && std::abs(d) < realsmall ){
        if(std::abs(c) < realsmall){
            // u.l. coincident zeros with l.r. of A=0; do nothing
            return 0;
        }
        else{
            // u.l. coincident zeros; put 0 in l.r. of A
            wz2 = arma::eye<arma::cx_mat>(2,2);
            //
            xy.set_size(1,2);
            xy(0,0) = c; xy(0,1) = - b;
            //
            xy = xy / arma::as_scalar(arma::sqrt(xy*arma::trans(xy)));
            xy1 = xy;
            xyt = arma::trans(xy);
            //
            xy2(0,0) = xyt(1,0);
            xy2(0,1) = -xyt(0,0);
            //
            xy2(1,arma::span()) = xy1;
        }
    }
    else{
        // usual case
        wz.set_size(1,2);
        wz(0,0) = c*e - f*b; wz(0,1) = std::conj(c*d - f*a);
        //
        xy.set_size(1,2);
        xy(0,0) = std::conj(b*d - e*a); xy(0,1) = std::conj(c*d - f*a);
        //
        n = arma::as_scalar(arma::sqrt(wz*arma::trans(wz)));
        m = arma::as_scalar(arma::sqrt(xy*arma::trans(xy)));
        //
        if( std::abs(m) < 1e-14){
            // all elements of A and B proportional
            return 0;
        }
        //
        wz = wz / n;
        xy = xy / m;
        //
        wz1 = wz; xy1 = xy;
        wzt = arma::trans(wz); xyt = arma::trans(xy);
        //
        wz2(0,arma::span()) = wz1;
        wz2(1,0) = - wzt(1,0); wz2(1,1) = wzt(0,0);
        //
        xy2(0,arma::span()) = xy1;
        xy2(1,0) = - xyt(1,0);
        xy2(1,1) = xyt(0,0);
        //
    }
    //
    wz = wz2;
    xy = xy2;
    //
    A(arma::span(i,i+1),arma::span()) = xy*A(arma::span(i,i+1),arma::span());
    B(arma::span(i,i+1),arma::span()) = xy*B(arma::span(i,i+1),arma::span());
    A(arma::span(),arma::span(i,i+1)) = A(arma::span(),arma::span(i,i+1))*wz;
    B(arma::span(),arma::span(i,i+1)) = B(arma::span(),arma::span(i,i+1))*wz;
    Z(arma::span(),arma::span(i,i+1)) = Z(arma::span(),arma::span(i,i+1))*wz;
    Q(arma::span(i,i+1),arma::span()) = xy*Q(arma::span(i,i+1),arma::span());
    //
    return 0;
    //
}
//
int qzdiv(double stake, arma::cx_mat &A, arma::cx_mat &B, arma::cx_mat &Q, arma::cx_mat &Z){
    //
    int n = A.n_rows;
    //
    arma::cx_colvec a = A.diag();
    arma::cx_colvec b = B.diag();
    //
    arma::mat roots(n,2);
    roots(arma::span(),0) = arma::abs(a);
    roots(arma::span(),1) = arma::abs(b);
    //
    int i,j,k,m;
    for(i = 0; i < n; ++i){
        if( roots(i,0) < 1.e-13 ){
            roots(i,0) += - (roots(i,0) + roots(i,1));
        }
    }
    //
    roots(arma::span(),1) = roots(arma::span(),1) / roots(arma::span(),0);
    //
    //
    //
    double tmp = 0;
    //
    for(i = n; i >= 1; --i){
        //
        m=0;
        //
        for(j = i; j >= 1; --j){
            if (roots(j-1,1) > stake || roots(j-1,1) < -0.1){
                m=j;
                break;
            }
        }
        if(m==0){
            return 0;
        }
        //
        for(k=m; k <= i-1; ++k){
            qzswitch(k-1,A,B,Q,Z); // notice that we have k - 1 here, not k as in Matlab code
            tmp = roots(k-1,1);
            roots(k-1,1) = roots(k,1);
            roots(k,1) = tmp;
        }
    }
    //
    return 0;
    
}
//
// Main gensys function, based on the discrete-time dynamic system:
//
//      Gamma0*y(t) = C + Gamma1*y(t-1) + Psi*z(t) + Pi*eta(t)
//
// where         z(t) is a vector of exogenous shocks,
//             eta(t) is a vector of 'one-step-ahead' expectational errors
//
SEXP gensysCpp(SEXP mGamma0, SEXP mGamma1, SEXP mC, SEXP mPsi, SEXP mPi)
{
#ifdef PREDEF_PLATFORM_UNIX
    try {
        //
        //Suppress warnings:
        std::ostream nullstream(0);
        arma::set_stream_err2(nullstream);
        //
        arma::mat Gamma0 = as<arma::mat>(mGamma0);
        arma::mat Gamma1 = as<arma::mat>(mGamma1);
        arma::mat C      = as<arma::mat>(mC);
        arma::mat Psi    = as<arma::mat>(mPsi);
        arma::mat Pi     = as<arma::mat>(mPi);
        //
        int n = Gamma0.n_rows;
        //
        arma::mat G1 = arma::zeros<arma::mat>(n,n);
        arma::mat Cons = arma::zeros<arma::mat>(n,1);
        arma::mat impact = arma::zeros<arma::mat>(n,Psi.n_cols);
        arma::cx_mat gev;
        arma::vec eu;
        //
        eu.set_size(2);
        eu.zeros();
        //
        arma::cx_mat xGamma0 = arma::zeros<arma::cx_mat>(n,n);
        arma::cx_mat xGamma1 = arma::zeros<arma::cx_mat>(n,n);
        //
        xGamma0.set_real(Gamma0);
        xGamma1.set_real(Gamma1);
        /*
        *
        * QZ decomposition: Gamma0 = Q**H S Z**H,  Gamma1 = Q**H T Z**H
        *
        */
        arma::cx_mat Q; arma::cx_mat Z; arma::cx_mat S; arma::cx_mat T;
        //
        arma::qz(S, T, Q, Z, xGamma0, xGamma1);
        //
        arma::cx_vec alpha_mat = S.diag();
        arma::cx_vec beta_mat = T.diag();
        //
        int nstable = 0; // this is 'nunstab' in Sims' code
        double zxz = 0;
        double divhat;
        double div = 1.01;
        double small = 1e-06;
        //
        int i;
        for(i = 0; i < n; ++i){
          //
          if( std::abs(alpha_mat(i)) > 0 ){
            divhat = std::abs(beta_mat(i))/std::abs(alpha_mat(i));
            //
            if( 1 + small < divhat && divhat <= div){
              div = 0.5*(1 + divhat);
            }
            //
          }
          //
          if( std::abs(beta_mat(i)) > div*std::abs(alpha_mat(i)) ){
            nstable = nstable + 1;
          }
          //
          if( std::abs(alpha_mat(i)) < small && std::abs(beta_mat(i)) < small ){
            zxz = 1;
          }
        }
        //
        if(zxz == 0){
          qzdiv(div,S,T,Q,Z);
        }
        else{ // zxz == 1
          eu(0) = -2; eu(1) = -2;
          //
          return Rcpp::List::create(Rcpp::Named("G1") = G1, Rcpp::Named("Cons") = Cons, Rcpp::Named("impact") = impact, Rcpp::Named("eu") = eu);
        }
        //
        gev.set_size(n,2);
        gev.col(0) = S.diag();
        gev.col(1) = T.diag();
        //
        arma::cx_mat Q1 = Q(arma::span(0,n-nstable-1),arma::span());
        arma::cx_mat Q2 = Q(arma::span(n-nstable,n-1),arma::span());
        arma::cx_mat Z1 = arma::trans(Z(arma::span(),arma::span(0,n-nstable-1)));
        arma::cx_mat Z2 = arma::trans(Z(arma::span(),arma::span(n-nstable,n-1)));
        arma::cx_mat S2 = S(arma::span(n-nstable,n-1),arma::span(n-nstable,n-1));
        arma::cx_mat T2 = T(arma::span(n-nstable,n-1),arma::span(n-nstable,n-1));
        //
        arma::cx_mat etawt = Q2*Pi;
        int neta = Pi.n_cols;
        //
        arma::cx_mat ueta = arma::zeros<arma::cx_mat>(0,0);
        arma::vec deta = arma::zeros<arma::vec>(0);
        arma::cx_mat veta = arma::zeros<arma::cx_mat>(neta,0);
        //
        if(nstable > 0){
          arma::cx_mat tueta; arma::vec tdeta; arma::cx_mat tveta;
          arma::uvec bigev;
          //
          arma::svd(tueta,tdeta,tveta,etawt);
          //
          bigev = arma::find(tdeta > small);
          //
          ueta = tueta.cols(bigev);
          veta = tveta.cols(bigev);
          deta = tdeta(bigev);
          //
          if(bigev.n_elem >= nstable){ // existence
            eu(0) = 1;
          }
        }
        //
        //
        //
        arma::cx_mat etawt1 = arma::zeros<arma::cx_mat>(0,neta);
        arma::cx_mat ueta1 = arma::zeros<arma::cx_mat>(0,0);
        arma::vec deta1 = arma::zeros<arma::vec>(0);
        arma::cx_mat veta1 = arma::zeros<arma::cx_mat>(neta,0);
        //
        if(nstable != n){
          arma::cx_mat tueta1; arma::vec tdeta1; arma::cx_mat tveta1;
          //
          etawt1 = Q1*Pi;
          //
          arma::svd(tueta1,tdeta1,tveta1,etawt1);
          //
          arma::uvec bigev2 = arma::find(tdeta1 > small);
          //
          ueta1 = tueta1.cols(bigev2);
          veta1 = tveta1.cols(bigev2);
          deta1 = tdeta1(bigev2);
        }
        //
        //
        //
        arma::cx_mat loose_temp;
        int uniq = 0;
        if(veta.n_rows==0){
          uniq = 1;
        }
        else{
          loose_temp = veta1 - veta*arma::trans(veta)*veta1;
          //
          arma::cx_mat ul; arma::vec dl; arma::cx_mat vl;
          arma::svd(ul,dl,vl,loose_temp);
          //
          arma::uvec kfind = arma::find( arma::abs(dl) > small*n);
          //
          if(kfind.n_elem == 0){
            uniq = 1;
          }
          else{
            uniq = 0;
          }
        }
        //
        if(uniq == 1){ // uniqueness
          eu(1) = 1;
        }
        /*
        *
        * Now put it all together
        *
        */
        arma::mat detamat = arma::diagmat(deta);   // put the singular values in diagonal matrix form
        arma::mat deta1mat = arma::diagmat(deta1);
        //
        arma::cx_mat tmat(n - nstable,n);
        tmat(arma::span(),arma::span(0,n-nstable-1)) = arma::eye<arma::cx_mat>(n-nstable,n-nstable);
        tmat(arma::span(),arma::span(n-nstable,n-1)) = - arma::trans((ueta*(arma::inv(detamat)*arma::trans(veta))*veta1*deta1mat*arma::trans(ueta1)));
        //
        //
        //
        arma::cx_mat G0(n,n);
        G0.zeros();
        G0(arma::span(0,n-nstable-1),arma::span()) = tmat * S;
        G0(arma::span(n-nstable,n-1),arma::span(n-nstable,n-1)) = arma::eye<arma::cx_mat>(nstable,nstable);
        //
        arma::cx_mat G0I = arma::inv(G0);
        //
        arma::cx_mat G1_temp(n,n);
        G1_temp.zeros();
        G1_temp(arma::span(0,n-nstable-1),arma::span()) = tmat * T;
        //
        G1_temp = G0I*G1_temp;
        //
        int usix = n - nstable + 1;
        //
        arma::cx_mat Cons_temp;
        bool Cstatus = arma::any(arma::vectorise(C));
        if(Cstatus==true){ // if any of the elements of 'C' are non-zero...
          arma::cx_mat C2(n,C.n_cols);
          C2.zeros();
          //
          C2(arma::span(0,n-nstable-1),arma::span()) = tmat * Q * C;
          C2(arma::span(n-nstable,n-1),arma::span()) = arma::inv( S(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) - T(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) ) * Q2 * C;
          //
          Cons_temp = G0I*C2;
        }else{ // ... otherwise set to zero. This avoids unnecessary calculations.
          Cons_temp = arma::zeros<arma::cx_mat>(n,C.n_cols);
        }
        //
        arma::cx_mat impact_temp(n,Psi.n_cols);
        impact_temp.zeros();
        impact_temp(arma::span(0,n-nstable-1),arma::span()) = tmat * Q * Psi;
        impact_temp = G0I * impact_temp;
        //
        //
        //
        G1 = arma::real(Z * G1_temp * arma::trans(Z));
        Cons = arma::real(Z * Cons_temp);
        impact = arma::real(Z * impact_temp);
        //
        return Rcpp::List::create(Rcpp::Named("G1") = G1,Rcpp::Named("Cons") = Cons,Rcpp::Named("impact") = impact,Rcpp::Named("eu") = eu);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: gensys C++ exception (unknown reason)" );
    }
    return R_NilValue;
#else
    ::Rf_error( "BMR: cannot use gensys without full LAPACK" );
    return R_NilValue;
#endif
}

int lyapunov_dbl (arma::mat &X, arma::mat F, int max_iter)
{
    //
    int j;
    //
    arma::mat X_O = X; // Old
    arma::mat X_N = X; // New
    //
    arma::mat ItMat = F;
    //
    for(j=1;j<=max_iter;j++){
        //
        X_N = X_O + ItMat*X_O*ItMat.t();
        //
        ItMat *= ItMat;
        X_O = X_N;
    }
    //
    X = X_N;
    //
    return 0;
}

SEXP DSGEKalman(SEXP dsgedata, SEXP mObserveMat, SEXP mObsCons, SEXP mF, SEXP mG, SEXP mshocks,
                 SEXP mR, SEXP mmax_iter)
{
    //
    try {
        //
        arma::mat data = as<arma::mat>(dsgedata);
        //
        arma::mat F = as<arma::mat>(mF);
        arma::mat G = as<arma::mat>(mG);
        arma::mat Q = as<arma::mat>(mshocks);
        //
        arma::mat C = as<arma::mat>(mObsCons);
        arma::mat H = as<arma::mat>(mObserveMat);
        arma::mat R = as<arma::mat>(mR);
        //
        int max_iter = as<int>(mmax_iter);
        //
        double loglikelihood = 0;
        //
        int n = data.n_rows;
        int k = data.n_cols;
        //
        double nlog2pi = k*log(2*arma::datum::pi);
        //
        int i;
        //
        arma::mat tH = H.t();
        arma::mat tF = F.t();
        //
        arma::mat GQG = G*Q*G.t();
        /*
         * Solve for the steady-state covariance matrix of the
         * state vector using a doubling algorithm
         */
        arma::mat InitialState(F.n_cols,1); InitialState.zeros();
        arma::mat InitialCov = GQG;
        //
        lyapunov_dbl(InitialCov, F, max_iter);
        /*
         * Now initialize the filter using the steady-state
         * covariance matrix
         */
        arma::mat StatePredicted = F*InitialState;
        arma::mat StateCovPredicted = F*InitialCov*tF + GQG;
        arma::mat Sigma = tH*StateCovPredicted*H + R;
        arma::mat invSigma = arma::inv_sympd(Sigma);
        //
        arma::mat KalmanGain = StateCovPredicted*H*invSigma;
        arma::mat KalmanResid = data.row(0).t() - tH*StatePredicted - C;
        arma::mat StateFiltered = StatePredicted + KalmanGain*KalmanResid;
        arma::mat StateCovFiltered = StateCovPredicted - KalmanGain*tH*StateCovPredicted;
        loglikelihood += nlog2pi + log(arma::det(Sigma)) + arma::as_scalar(KalmanResid.t()*invSigma*KalmanResid);
        //
        for(i=2; i<=n; i++){
            StatePredicted = F*StateFiltered;
            StateCovPredicted = F*StateCovFiltered*tF + GQG;
            Sigma = tH*StateCovPredicted*H + R;
            invSigma = arma::inv_sympd(Sigma);
            //
            KalmanGain = StateCovPredicted*H*invSigma;
            KalmanResid = data.row(i-1).t() - tH*StatePredicted - C;
            StateFiltered = StatePredicted + KalmanGain*KalmanResid;
            StateCovFiltered = StateCovPredicted - KalmanGain*tH*StateCovPredicted;
            loglikelihood += nlog2pi + log(arma::det(Sigma)) + arma::as_scalar(KalmanResid.t()*invSigma*KalmanResid);
        }
        // Returns the negative of the loglikelihood
        loglikelihood = loglikelihood/2;
        //
        //
        return Rcpp::List::create(Rcpp::Named("dsgelike") = loglikelihood);
    } catch( std::exception &ex ) {
        //forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: DSGE Kalman C++ exception (unknown reason)" );
    }
    //return R_NilValue;
    double LLKerr = 1000000;
    return Rcpp::List::create(Rcpp::Named("dsgelike") = LLKerr);
}

// This returns the period-by-period likelihood value and filtered values of the state variables.
SEXP DSGEKalmanFilt(SEXP dsgedata, SEXP mObserveMat,
                    SEXP mObsCons, SEXP mF, SEXP mG, SEXP mshocks,
                    SEXP mR, SEXP mmax_iter)
{
    try {
        //
        arma::mat data = as<arma::mat>(dsgedata);
        //
        arma::mat F = as<arma::mat>(mF);
        arma::mat G = as<arma::mat>(mG);
        arma::mat Q = as<arma::mat>(mshocks);
        //
        arma::mat C = as<arma::mat>(mObsCons);
        arma::mat H = as<arma::mat>(mObserveMat);
        arma::mat R = as<arma::mat>(mR);
        //
        int max_iter = as<int>(mmax_iter);
        //
        //
        int n = data.n_rows;
        int k = data.n_cols;
        //
        double nlog2pi = k*log(2*arma::datum::pi);
        //
        int i;
        //
        arma::mat tH = H.t();
        arma::mat tF = F.t();
        //
        arma::mat loglikelihood(n,1); loglikelihood.zeros();
        arma::mat StateMat = arma::zeros<arma::mat>(F.n_cols,n);
        //
        arma::mat GQG = G*Q*G.t();
        /*
         * Solve for the steady-state covariance matrix of the
         * state vector using a doubling algorithm
         */
        arma::mat InitialState(F.n_cols,1); InitialState.zeros();
        arma::mat InitialCov = GQG;
        //
        lyapunov_dbl(InitialCov, F, max_iter);
        /*
         * Now initialize the filter using the steady-state
         * covariance matrix
         */
        arma::mat StatePredicted = F*InitialState;
        arma::mat StateCovPredicted = F*InitialCov*tF + GQG;
        arma::mat Sigma = tH*StateCovPredicted*H + R;
        arma::mat invSigma = arma::inv_sympd(Sigma);
        //
        arma::mat KalmanGain = StateCovPredicted*H*invSigma;
        arma::mat KalmanResid = data.row(0).t() - tH*StatePredicted - C;
        arma::mat StateFiltered = StatePredicted + KalmanGain*KalmanResid;
        StateMat.col(0) = StateFiltered;
        arma::mat StateCovFiltered = StateCovPredicted - KalmanGain*tH*StateCovPredicted;
        loglikelihood.row(0) = nlog2pi + log(arma::det(Sigma)) + KalmanResid.t()*invSigma*KalmanResid;
        //
        for(i=2; i<=n; i++){
            StatePredicted = F*StateFiltered;
            StateCovPredicted = F*StateCovFiltered*tF + GQG;
            Sigma = tH*StateCovPredicted*H + R;
            invSigma = arma::inv_sympd(Sigma);
            //
            KalmanGain = StateCovPredicted*H*invSigma;
            KalmanResid = data.row(i-1).t() - tH*StatePredicted - C;
            StateFiltered = StatePredicted + KalmanGain*KalmanResid;
            StateMat.col(i-1) = StateFiltered;
            StateCovFiltered = StateCovPredicted - KalmanGain*tH*StateCovPredicted;
            loglikelihood.row(i-1) = nlog2pi + log(arma::det(Sigma)) + KalmanResid.t()*invSigma*KalmanResid;
        }
        // Returns the negative of the loglikelihood
        loglikelihood = 0.5*loglikelihood;
        //
        return Rcpp::List::create(Rcpp::Named("dsgelike") = loglikelihood,Rcpp::Named("dsgestate") = StateMat);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: DSGE Kalman Filter C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP DSGECR(SEXP dsgedata, SEXP mObserveMat,
            SEXP mObsCons, SEXP mF, SEXP mG, SEXP mshocks,
            SEXP mR, SEXP mmax_iter)
{
    //
    try {
        //
        arma::mat data = as<arma::mat>(dsgedata);
        //
        arma::mat F = as<arma::mat>(mF);
        arma::mat G = as<arma::mat>(mG);
        arma::mat Q = as<arma::mat>(mshocks);
        //
        arma::mat C = as<arma::mat>(mObsCons);
        arma::mat H = as<arma::mat>(mObserveMat);
        arma::mat R = as<arma::mat>(mR);
        //
        int max_iter = as<int>(mmax_iter);
        //
        double loglikelihood = 0;
        //
        int n = data.n_rows;
        int k = data.n_cols;
        //
        double nlog2pi = k*log(2*arma::datum::pi);
        //
        int i;
        //
        arma::mat tH = H.t();
        arma::mat tF = F.t();
        //
        arma::mat GQG = G*Q*G.t();
        /*
         * Solve for the steady-state covariance matrix of the
         * state vector using a doubling algorithm
         */
        arma::mat InitialState(F.n_cols,1); InitialState.zeros();
        arma::mat InitialCov = GQG;
        //
        lyapunov_dbl(InitialCov, F, max_iter);
        /*
         * Now initialize the filter using the steady-state
         * covariance matrix
         */
        arma::mat StateFiltered = InitialState;
        arma::mat StateCovPredicted = F*InitialCov*tF + GQG;
        arma::mat Sigma = tH*StateCovPredicted*H + R;
        arma::mat iSigma = arma::inv_sympd(Sigma);
        //
        arma::mat St = F*StateCovPredicted*H;
        arma::mat Mt = -iSigma;
        arma::mat Kt = St*iSigma;
        //
        arma::mat Resid = data.row(0).t() - tH*StateFiltered - C;
        //
        loglikelihood += nlog2pi + log(arma::det(Sigma)) + arma::as_scalar(Resid.t()*iSigma*Resid);
        //
        StateFiltered = F*StateFiltered + Kt*Resid;
        //
        arma::mat tHSt = tH*St;
        arma::mat MSpZp = Mt*trans(tHSt);
        arma::mat FSt = F*St;
        //
        arma::mat Sigma1  = Sigma + tHSt*MSpZp;
        arma::mat iSigma1 = arma::inv_sympd(Sigma1);
        //
        Kt = (Kt*Sigma + FSt*MSpZp)*iSigma1;
        St = FSt - Kt*tHSt;
        Mt += MSpZp*iSigma*trans(MSpZp);
        Sigma = Sigma1;
        iSigma = iSigma1;
        //
        for(i=2; i<=n; i++){
            //
            Resid = data.row(i-1).t() - tH*StateFiltered - C;
            //
            loglikelihood += nlog2pi + log(arma::det(Sigma)) + arma::as_scalar(Resid.t()*iSigma*Resid);
            //
            StateFiltered = F*StateFiltered + Kt*Resid;
            //
            tHSt = tH*St;
            MSpZp = Mt*trans(tHSt);
            FSt = F*St;
            //
            Sigma1 += tHSt*MSpZp;
            iSigma1 = arma::inv_sympd(Sigma1);
            //
            Kt = (Kt*Sigma + FSt*MSpZp)*iSigma1;
            St = FSt - Kt*tHSt;
            Mt += MSpZp*iSigma*trans(MSpZp);
            Sigma = Sigma1;
            iSigma = iSigma1;
        }
        // Returns the negative of the loglikelihood
        loglikelihood = loglikelihood/2;
        //
        return Rcpp::List::create(Rcpp::Named("dsgelike") = loglikelihood);
    } catch( std::exception &ex ) {
        //  forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: DSGE CR C++ exception (unknown reason)" );
    }
    //return R_NilValue;
    double LLKerr = 1000000;
    return Rcpp::List::create(Rcpp::Named("dsgelike") = LLKerr);
}