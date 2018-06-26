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
 * simulate from a VAR model: Y = X * beta + e
 */

inline
arma::mat
var_sim(const arma::mat& var_coefs, const bool cons_term, const int n_out, const int n_burnin)
{
    const int M = var_coefs.n_cols;
    const int K = var_coefs.n_rows;

    const int c_int = !!cons_term;
    const int p = (K - c_int) / M;

    //
    
    arma::mat Y = arma::zeros(n_out + p + n_burnin,M);

    arma::mat mean_vec = arma::zeros(M,1);

    if (cons_term) {
        arma::mat Phi = var_coefs.row(0);

        arma::mat poly_mat = arma::eye(M,M);
        
        for (int i=0; i < p; i++) {
            poly_mat -= var_coefs.rows(1 + i*M,(i+1)*M);
        }

        mean_vec = Phi*arma::inv(poly_mat);
    }

    //

    Y.rows(0,p) = arma::repmat(mean_vec,p+1,1);

    for (int i=p; i < (n_out + p + n_burnin); i++) {
        
        if (cons_term) {
            Y.row(i) = var_coefs.row(0);
        }

        for (int j=0; j < p; j++) {
            Y.row(i) += Y.row(i-j-1) * var_coefs.rows(c_int + j*M, (j+1)*M);
        }

        Y.row(i) += arma::trans(stats::rmvnorm<arma::mat,double>(arma::zeros(M,1),arma::eye(M,M),true));
    }
    
    Y.shed_rows(0,n_burnin-1);

    return Y;
}
