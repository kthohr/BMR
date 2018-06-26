/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the BM++ C++ library.
  ##
  ##   BM++ is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   BM++ is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with BM++. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/* 
 *  Main gensys function, based on a discrete-time dynamic system of the form
 *
 *         Gamma0*y(t) = GammaC + Gamma1*y(t-1) + Psi*z(t) + Pi*eta(t)
 *
 *  where:         z(t) is a vector of exogenous shocks,
 *               eta(t) is a vector of 'one-step-ahead' expectational errors
 *
 *  Output is of the form:
 *
 *         y(t) = Cons + G1*y(t-1) + impact*z(t) + ywt*inv(I - fmat*inv(L))*fwt*z(t+1).
 *
 *  In a lot of cases z is iid, and so the final term drops out.
 *
 * Original Matlab code by Chris Sims
 *
 * Ported and adapted to C++ by
 * Keith O'Hara
 * 06/29/15
 *
 * This version:
 * 08/16/17
 */

#include "misc/misc.hpp"
namespace bm {
    #include "lres/gensys_solver.hpp"
}

int bm::gensys_solver(const arma::mat& Gamma0, const arma::mat& Gamma1, const arma::mat& GammaC, const arma::mat& Psi, const arma::mat& Pi,
                      arma::mat& G1, arma::mat& Cons, arma::mat& impact, 
                      arma::cx_mat* fmat_out, arma::cx_mat* fwt_out, arma::cx_mat* ywt_out, arma::cx_mat* gev_out, arma::vec* eu_out, arma::mat* loose_out)
{
    //
    const double small_val = 1e-06;

    const int n = Gamma0.n_rows;
    
    if (eu_out) {
        eu_out->zeros(2);
    }
    
    arma::cx_mat xGamma0 = arma::zeros<arma::cx_mat>(n,n);
    arma::cx_mat xGamma1 = arma::zeros<arma::cx_mat>(n,n);
    
    xGamma0.set_real(Gamma0);
    xGamma1.set_real(Gamma1);
    
    //
    // QZ decomposition: Gamma0 = Q**H S Z**H,  Gamma1 = Q**H T Z**H

    arma::cx_mat Q; arma::cx_mat Z; arma::cx_mat S; arma::cx_mat T;
    
    arma::qz(S, T, Q, Z, xGamma0, xGamma1);
    
    arma::cx_vec alpha_mat = S.diag();
    arma::cx_vec beta_mat = T.diag();

    //
    
    int n_stable = 0; // this is 'nunstab' in Sims' code

    double zxz = 0;
    double divhat;
    double div = 1.01;
    
    for (int i = 0; i < n; i++) {
        //
        if (std::abs(alpha_mat(i)) > 0) {
            divhat = std::abs(beta_mat(i))/std::abs(alpha_mat(i));
            
            if (1 + small_val < divhat && divhat <= div) {
                div = 0.5*(1 + divhat);
            }
        }
        
        if (std::abs(beta_mat(i)) > div*std::abs(alpha_mat(i))) {
            n_stable += 1;
        }
        
        if (std::abs(alpha_mat(i)) < small_val && std::abs(beta_mat(i)) < small_val) {
            zxz = 1;
        }
    }

    // error check
    if (zxz == 0) {
        qz_div(div,S,T,Q,Z);
    } else { // zxz == 1
        std::cout << "Coincident zeros. Indeterminacy and/or nonexistence." << std::endl;

        if (eu_out) {
            (*eu_out)(0) = -2;
            (*eu_out)(1) = -2;
        }

        // set solution objects to zero (avoids segfault issues in estimation)
        G1.zeros(n,n);
        Cons.zeros(n,1);
        impact.zeros(n,Psi.n_cols);
        //
        return -2;
    }

    if (gev_out) {
        gev_out->set_size(n,2);
        gev_out->col(0) = S.diag();
        gev_out->col(1) = T.diag();
    }

    //
    // Setup Q and Z matrices

    arma::cx_mat Q1, Q2, Z1, Z2, S2, T2;
    
    if (n_stable < n) {
        Q1 = Q.rows(0,n-n_stable-1);
        Z1 = arma::trans(Z.cols(0,n-n_stable-1));
        
        if (n_stable != 0) {
            Q2 = Q.rows(n-n_stable,n-1);
            Z2 = arma::trans(Z.cols(n-n_stable,n-1));
            
            S2 = S(arma::span(n-n_stable,n-1),arma::span(n-n_stable,n-1));
            T2 = T(arma::span(n-n_stable,n-1),arma::span(n-n_stable,n-1));
        } else {
            Q2.set_size(0,n);
            Z2.set_size(0,n);
            
            S2.set_size(0,0);
            T2.set_size(0,0);
        }
    } else {
        Q1.set_size(0,n);
        Q2 = Q.rows(n-n_stable,n-1);
        
        Z1.set_size(0,n);
        Z2 = arma::trans(Z.cols(n-n_stable,n-1));
        
        S2 = S(arma::span(n-n_stable,n-1),arma::span(n-n_stable,n-1));
        T2 = T(arma::span(n-n_stable,n-1),arma::span(n-n_stable,n-1));
    }
    //
    arma::cx_mat etawt = Q2*Pi;
    int n_eta = Pi.n_cols;
    
    arma::cx_mat ueta = arma::zeros<arma::cx_mat>(0,0);
    arma::vec deta = arma::zeros<arma::vec>(0);
    arma::cx_mat veta = arma::zeros<arma::cx_mat>(n_eta,0);
    //
    if (n_stable > 0) {
        arma::cx_mat tueta; arma::vec tdeta; arma::cx_mat tveta;
        arma::uvec bigev;
        
        arma::svd(tueta,tdeta,tveta,etawt);
        //
        bigev = arma::find(tdeta > small_val);
        
        ueta = tueta.cols(bigev);
        veta = tveta.cols(bigev);
        deta = tdeta(bigev);
        //
        if (eu_out && (int) bigev.n_elem >= n_stable) { // existence
            (*eu_out)(0) = 1;
        }
    }
    
    //
    
    arma::cx_mat etawt1 = arma::zeros<arma::cx_mat>(0,n_eta), ueta1 = arma::zeros<arma::cx_mat>(0,0), veta1 = arma::zeros<arma::cx_mat>(n_eta,0);
    arma::vec deta1 = arma::zeros<arma::vec>(0);
    
    if (n_stable != n) {
        arma::cx_mat tueta1; arma::vec tdeta1; arma::cx_mat tveta1;
        //
        etawt1 = Q1*Pi;
        //
        arma::svd(tueta1,tdeta1,tveta1,etawt1);
        //
        arma::uvec bigev2 = arma::find(tdeta1 > small_val);
        //
        ueta1 = tueta1.cols(bigev2);
        veta1 = tveta1.cols(bigev2);
        deta1 = tdeta1(bigev2);
    }
    
    //
    // uniqueness check
    
    int uniq = 0, n_loose = 0;
    // arma::cx_mat loose_temp;

    if (veta.n_rows == 0) {
        uniq = 1;
    } else {
        arma::cx_mat loose_temp = veta1 - veta*arma::trans(veta)*veta1;
        //
        arma::cx_mat ul; arma::vec dl; arma::cx_mat vl;
        arma::svd(ul,dl,vl,loose_temp);
        //
        arma::uvec kfind = arma::find( arma::abs(dl) > small_val*n);
        //
        if (kfind.n_elem == 0) {
            uniq = 1;
        } else {
            uniq = 0;
        }

        n_loose = kfind.n_elem;
    }
    
    if (uniq == 1) { // uniqueness
        if (eu_out) {
            (*eu_out)(1) = 1;
        }
    } else {
        std::cout << "Indeterminacy. " << n_loose << " loose endog errors." << std::endl;
    }
    
    //
    // Now put it all together
    
    arma::mat deta_mat  = arma::diagmat(deta);   // put the singular values in diagonal matrix form
    arma::mat deta1_mat = arma::diagmat(deta1);
    
    arma::cx_mat tmat(n - n_stable,n);

    tmat.cols(0,n-n_stable-1) = arma::eye<arma::cx_mat>(n-n_stable,n-n_stable);

    if (n_stable != 0) {
        tmat.cols(n-n_stable,n-1) = - arma::trans( ueta*(arma::inv(deta_mat)*arma::trans(veta)) * veta1*deta1_mat*arma::trans(ueta1) );
    }
    
    arma::cx_mat G0(n,n);

    G0.zeros();
    G0.rows(0,n-n_stable-1) = tmat * S;

    if (n_stable != 0) {
        G0(arma::span(n-n_stable,n-1),arma::span(n-n_stable,n-1)) = arma::eye<arma::cx_mat>(n_stable,n_stable);
    }
    //
    arma::cx_mat G0_inv = arma::inv(G0);
    
    arma::cx_mat G1_temp(n,n);
    G1_temp.zeros();
    G1_temp.rows(0,n-n_stable-1) = tmat * T;
    
    G1_temp = G0_inv*G1_temp;
    
    //
    // intercept term

    int usix = n - n_stable + 1;
    
    arma::cx_mat Cons_temp;
    bool Cstatus = arma::any(arma::vectorise(GammaC)); // if any of the elements of 'GammaC' are non-zero...

    if ( Cstatus ) {
        arma::cx_mat C2(n,GammaC.n_cols);
        C2.zeros();
        //
        C2.rows(0,n-n_stable-1) = tmat * Q * GammaC;

        if (n_stable != 0) {
            C2.rows(n-n_stable,n-1) = arma::inv( S(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) - T(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) ) * Q2 * GammaC;
        }
        //
        Cons_temp = G0_inv*C2;
    } else { // ... otherwise set to zero. This avoids unnecessary calculations.
        Cons_temp = arma::zeros<arma::cx_mat>(n,GammaC.n_cols);
    }

    //
    // impack matrix

    arma::cx_mat impact_temp(n,Psi.n_cols);
    impact_temp.zeros();

    impact_temp.rows(0,n-n_stable-1) = tmat * Q * Psi;
    impact_temp = G0_inv * impact_temp;
    
    //
    // additional output

    if (n_stable != 0) {
        if (fmat_out) {
            *fmat_out = arma::inv( T(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) ) * S(arma::span(usix-1,n-1),arma::span(usix-1,n-1));
        }

        if (fwt_out) {
            *fwt_out = - arma::inv( T(arma::span(usix-1,n-1),arma::span(usix-1,n-1)) ) * Q2 * Psi;
        }
        //
        if (ywt_out) {
            *ywt_out = Z * G0_inv.cols(usix-1,n-1);
        }
    } 
    // else {
    //     fmat.set_size(0,0);
    //     fwt.set_size(0,0);
    //     ywt.set_size(0,0);
    // }
    //

    if (loose_out) {
        // loose_temp.set_size(n,n_eta);
        arma::cx_mat loose_temp(n,n_eta);
        loose_temp.zeros();

        loose_temp.rows(0,n-n_stable-1) = etawt1 * (arma::eye<arma::cx_mat>(n_eta,n_eta) - veta * arma::trans(veta));
        loose_temp = G0_inv*loose_temp;

        *loose_out = arma::real(Z * loose_temp);
    }
    
    //
    // finish solution matrices

    G1 = arma::real(Z * G1_temp * arma::trans(Z));
    Cons = arma::real(Z * Cons_temp);
    impact = arma::real(Z * impact_temp);
    
    //
    return 0;
}

//
// shortened solver; only return main solver elements.

int bm::gensys_solver(const arma::mat& Gamma0, const arma::mat& Gamma1, const arma::mat& GammaC, const arma::mat& Psi, const arma::mat& Pi,
                       arma::mat& G1, arma::mat& Cons, arma::mat& impact)
{
    return gensys_solver(Gamma0,Gamma1,GammaC,Psi,Pi,G1,Cons,impact,nullptr,nullptr,nullptr,nullptr,nullptr,nullptr);
}

//
// internal functions

void bm::qz_switch(int i, arma::cx_mat& A, arma::cx_mat& B, arma::cx_mat& Q, arma::cx_mat& Z)
{
    //
    const double realsmall_val = 1e-08;
    //
    std::complex<double> a = A(i,i);         std::complex<double> d = B(i,i);
    std::complex<double> b = A(i,i+1);     std::complex<double> e = B(i,i+1);
    std::complex<double> c = A(i+1,i+1); std::complex<double> f = B(i+1,i+1);
    
    arma::cx_mat wz(2,1); arma::cx_mat wz1(2,1); arma::cx_mat wz2(2,2); arma::cx_mat wzt;
    arma::cx_mat xy(1,2); arma::cx_mat xy1; arma::cx_mat xy2(2,2); arma::cx_mat xyt;

    //

    std::complex<double> n,m;
    
    if (std::abs(c) < realsmall_val && std::abs(f) < realsmall_val) {
        if (std::abs(a) < realsmall_val) {
            // l.r. coincident 0's with u.l. of A=0; do nothing
            return;
        }
        else{
            // l.r. coincident zeros; put 0 in u.l. of a
            wz.set_size(2,1);
            
            wz(0,0) = b;
            wz(1,0) = -a;
            //
            wz = wz / arma::as_scalar(arma::sqrt(arma::trans(wz)*wz));
            wz1 = wz;
            wzt = arma::trans(wz);
            //
            wz2.col(0) = wz1;
            //
            wz2(0,1) = wzt(0,1);
            wz2(1,1) = -wzt(0,0);
            //
            xy2 = arma::eye<arma::cx_mat>(2,2);
        }
    }
    else if (std::abs(a) < realsmall_val && std::abs(d) < realsmall_val) {
        if (std::abs(c) < realsmall_val) {
            // u.l. coincident zeros with l.r. of A=0; do nothing
            return;
        } else {
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
            xy2.row(1) = xy1;
        }
    } else {
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
        if (std::abs(m) < 1e-14) {
            // all elements of A and B proportional
            return;
        }
        //
        wz = wz / n;
        xy = xy / m;
        //
        wz1 = wz; xy1 = xy;
        wzt = arma::trans(wz); xyt = arma::trans(xy);
        //
        wz2.row(0) = wz1;
        wz2(1,0) = - wzt(1,0); wz2(1,1) = wzt(0,0);
        //
        xy2.row(0) = xy1;
        xy2(1,0) = - xyt(1,0);
        xy2(1,1) = xyt(0,0);
        //
    }
    //
    wz = wz2;
    xy = xy2;
    //
    A.rows(i,i+1) = xy*A.rows(i,i+1);
    B.rows(i,i+1) = xy*B.rows(i,i+1);
    A.cols(i,i+1) = A.cols(i,i+1)*wz;
    B.cols(i,i+1) = B.cols(i,i+1)*wz;
    Z.cols(i,i+1) = Z.cols(i,i+1)*wz;
    Q.rows(i,i+1) = xy*Q.rows(i,i+1);
    //
}

void bm::qz_div(double stake, arma::cx_mat& A, arma::cx_mat& B, arma::cx_mat& Q, arma::cx_mat& Z)
{
    //
    int n = A.n_rows;
    //
    arma::cx_colvec a = A.diag();
    arma::cx_colvec b = B.diag();
    //
    arma::mat roots(n,2);
    roots.col(0) = arma::abs(a);
    roots.col(1) = arma::abs(b);
    //
    int i,j,k,m;
    for(i = 0; i < n; ++i){
        if( roots(i,0) < 1.e-13 ){
            roots(i,0) += - (roots(i,0) + roots(i,1));
        }
    }
    //
    roots.col(1) = roots.col(1) / roots.col(0);
    //
    double tmp = 0;
    
    for(i = n; i >= 1; --i){
        //
        m = 0;
        
        for (j = i; j >= 1; --j) {
            if (roots(j-1,1) > stake || roots(j-1,1) < -0.1) {
                m = j;
                break;
            }
        }

        if (m == 0) {
            return;
        }
        //
        for (k=m; k <= i-1; ++k) {
            qz_switch(k-1,A,B,Q,Z); // we have k - 1 here; not k as in Matlab code

            tmp = roots(k-1,1);
            roots(k-1,1) = roots(k,1);
            roots(k,1) = tmp;
        }
    }
    //
}
