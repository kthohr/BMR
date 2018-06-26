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
 * Uhlig's method for solving linear rational  expectations
 * models using a generalized eigen decomposition.
 */

#include "misc/misc.hpp"
namespace bm {
    #include "lres/uhlig_solver.hpp"
}

int bm::uhlig_solver(const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& D,
                     const arma::mat& F, const arma::mat& G, const arma::mat& H, const arma::mat& J, const arma::mat& K, const arma::mat& L, const arma::mat& M, const arma::mat& N,
                     const arma::vec* which_eig, arma::mat& P, arma::mat& Q, arma::mat& R, arma::mat& S, arma::cx_vec* eigen_vals_out, arma::cx_mat* eigen_vecs_out)
{
    const double bignum = 1e+08;

    //
    
    const int l = C.n_rows;
    const int n = C.n_cols;
    const int k = std::min(N.n_rows,N.n_cols);
    
    const int m = (l == 0) ? F.n_cols : A.n_cols;

    const bool pick_eig = (which_eig) ? ( (which_eig->n_elem > 0) ? true : false ) : false;
    
    //
    // setup

    arma::mat Xi(2*m,2*m);
    arma::mat Delta(2*m,2*m);
    
    arma::mat Psi = F;
    arma::mat Gamma = -G;
    arma::mat Theta = -H;
    
    arma::mat C_pinv = C;

    if (l > 0) {
        C_pinv = arma::pinv(C);
        Psi = F - J*C_pinv*A;
        Gamma = J*C_pinv*B - G + K*C_pinv*A;
        Theta = K*C_pinv*B - H;
    }
    
    Xi(arma::span(0,m-1),arma::span(0,m-1)) = Gamma;
    Xi(arma::span(0,m-1),arma::span(m,2*m-1)) = Theta;
    Xi(arma::span(m,2*m-1),arma::span(0,m-1)) = arma::eye(m,m);
    Xi(arma::span(m,2*m-1),arma::span(m,2*m-1)).zeros();
    
    Delta(arma::span(0,m-1),arma::span(0,m-1)) = Psi;
    Delta(arma::span(0,m-1),arma::span(m,2*m-1)).zeros();
    Delta(arma::span(m,2*m-1),arma::span(0,m-1)).zeros();
    Delta(arma::span(m,2*m-1),arma::span(m,2*m-1)) = arma::eye(m,m);

    //
    // Next perform a generalized eigen decomp:
    
    arma::cx_vec EigValue;
    arma::cx_mat EigVec;
    
    arma::eig_pair(EigValue, EigVec, Xi, Delta);
    
    arma::vec EigValueAbs = arma::abs(EigValue);
    arma::vec EigValueReal = arma::real(EigValue);

    //
    // Deal with infinite values (otherwise sort will return an error):
    
    // arma::uvec pos_inf_ind = arma::find_nonfinite(EigValueAbs);
    arma::uvec pos_inf_ind = arma::find(EigValueReal > bignum);
    arma::uvec neg_inf_ind = arma::find(EigValueReal < -bignum);

    int infin2 = pos_inf_ind.n_elem;
    int infin3 = neg_inf_ind.n_elem;
    
    if (infin2 > 0) {
        arma::cx_vec BigNum(infin2);
        BigNum.ones();
        
        EigValueAbs.elem(pos_inf_ind).fill( bignum );
        EigValue.elem(pos_inf_ind) = BigNum*bignum;
    }

    if (infin3 > 0) {
        arma::cx_vec BigNum(infin3);
        BigNum.ones();
        
        EigValue.elem(neg_inf_ind) = - BigNum*bignum;
    }

    //
    // Now sort the eigenvalues and eigenvectors, and calculate P
    
    arma::uvec indices = arma::sort_index(EigValueAbs);

    // ... from smallest to largest in absolute value
    arma::vec EigValueAbsSorted = EigValueAbs.elem(indices);
    arma::cx_vec EigValueSorted = EigValue.elem(indices);
    arma::cx_mat EigVecSorted = EigVec.cols(indices);
    
    // If the user prefers to 'choose' which eigenvalues to use...

    if (pick_eig) {
        // by 'egvecind-1' this implies that the elements of 'which_eig' begin at 1 and not 0
        indices = arma::conv_to<arma::uvec>::from(*which_eig - 1);
        
        EigValueSorted = EigValueSorted.elem(indices);
        EigVecSorted = EigVecSorted.cols(indices);
    }

    //

    if (eigen_vals_out) {
        *eigen_vals_out = EigValueSorted;
    }

    if (eigen_vecs_out) {
        *eigen_vecs_out = EigVecSorted;
    }
    
    arma::cx_vec LambdaVec = EigValueSorted.rows(0,m-1);
    arma::cx_mat Lambda = arma::diagmat(LambdaVec);
    
    arma::cx_mat Omega = EigVecSorted(arma::span(m,2*m-1),arma::span(0,m-1));

    P = arma::real( Omega*Lambda*arma::pinv(Omega) );

    //
    // Now calculate Q, R, and S

    C_pinv = C;
    
    arma::colvec LNMStacked = arma::vectorise(L*N + M);
    
    arma::mat V = arma::kron(N.t(),F) + arma::kron(arma::eye(k,k),(F*P + G));

    arma::vec QS;
    R.set_size(0,P.n_cols); R.zeros();

    if (l == 0) {
        QS = -arma::inv(V)*LNMStacked;
    } else {

        //
        // First R

        C_pinv = arma::pinv(C);
        R = - C_pinv*(A*P + B);
        
        //
        // then Q
        
        arma::mat V2 = arma::zeros((k*A.n_rows) + V.n_rows, (k*A.n_rows) + V.n_rows);
        
        V2(arma::span(0,(k*A.n_rows)-1),arma::span(0,(k*A.n_cols)-1)) = arma::kron(arma::eye(k,k),A);
        V2(arma::span(0,(k*A.n_rows)-1),arma::span(k*A.n_cols,(k*A.n_cols)+(k*C.n_cols)-1)) = arma::kron(arma::eye(k,k),C);
        V2(arma::span(k*A.n_rows,(k*A.n_rows)+V.n_rows-1),arma::span(0,(k*F.n_cols)-1)) = arma::kron(N.t(),F) + arma::kron(arma::eye(k,k),(F*P + J*R + G));
        V2(arma::span(k*A.n_rows,(k*A.n_rows)+V.n_rows-1),arma::span(k*F.n_cols,(k*F.n_cols)+(k*J.n_cols)-1)) = arma::kron(N.t(),J) + arma::kron(arma::eye(k,k),K);
        
        arma::vec DLNM = arma::join_cols(arma::vectorise(D),LNMStacked);

        //
        QS = - arma::inv(V2)*DLNM;
    }

    Q = arma::reshape(QS.rows(0,(m*k)-1), m, k);

    //
    // Finally, S...
    
    S.set_size(0,Q.n_cols); S.zeros();

    if (l > 0) {
        S = arma::reshape(QS.rows(m*k,(m+n)*k-1), n, k);
    }

    //

    return 0;
}

// without which_eig or returning the eigenvalues and eigenvectors
int bm::uhlig_solver(const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& D,
                     const arma::mat& F, const arma::mat& G, const arma::mat& H, const arma::mat& J, const arma::mat& K, const arma::mat& L, const arma::mat& M, const arma::mat& N,
                     arma::mat& P, arma::mat& Q, arma::mat& R, arma::mat& S)
{
    return uhlig_solver(A,B,C,D,F,G,H,J,K,L,M,N,nullptr,P,Q,R,S,nullptr,nullptr);
}
