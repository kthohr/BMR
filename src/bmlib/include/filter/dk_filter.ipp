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
 * Durbin-Koopman simulation smoother
 */

inline
arma::mat
dk_filter(const arma::mat& Y, const arma::mat& Z, const arma::mat& Q_draw, const arma::mat& Q_chol, const arma::mat& Sigma_draw, const arma::mat& Sigma_chol, const int n_adj, const int M, const int K)
{

    arma::mat eye_K = arma::eye(K*M,K*M);

    arma::mat beta_plus = arma::zeros(K*M,n_adj+1); // ensures that beta_1 ~ N(0,Q)
    arma::mat Y_plus = arma::zeros(M,n_adj);

    for (int i=0; i < n_adj; i++) {
        arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

        Y_plus.col(i)  = Z_t*beta_plus.col(i) + Sigma_chol*arma::randn(M,1);
        beta_plus.col(i+1) = beta_plus.col(i) + Q_chol*arma::randn(K*M,1);
    }

    arma::mat Y_star = Y - Y_plus;

    // Kalman filter

    arma::mat Lt_mat = arma::zeros(K*M*n_adj,K*M);
    arma::mat Ft_mat = arma::zeros(M*n_adj,M);

    arma::mat beta_filt = arma::zeros(K*M,n_adj+1);
    arma::mat kal_resid = arma::zeros(M,n_adj);
    arma::mat state_cov_filt = arma::zeros(K*M,K*M);

    for (int i=0; i < n_adj; i++) {
        arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

        //

        arma::mat Ft = Z_t*state_cov_filt*Z_t.t() + Sigma_draw;
        arma::mat inv_Ft = arma::inv_sympd(Ft);

        arma::mat kal_gain = state_cov_filt*Z_t.t()*inv_Ft;
        kal_resid.col(i) = Y_star.col(i) - Z_t*beta_filt.col(i);

        arma::mat Lt = eye_K - kal_gain*Z_t;

        beta_filt.col(i+1) = beta_filt.col(i) + kal_gain*kal_resid.col(i);
        state_cov_filt = state_cov_filt*Lt.t() + Q_draw;

        // storage

        Ft_mat.rows(i*M,(i+1)*M-1) = inv_Ft;
        Lt_mat.rows(i*K*M,(i+1)*K*M-1) = Lt;
    }

    // smoothed disturbance; eq (5) in Durbin and Koopman

    arma::mat rt_mat = arma::zeros(K*M,n_adj+1);

    for (int j=n_adj; j >= 1; j--) {
        arma::mat Z_t = Z.rows((j-1)*M,(j*M)-1);

        // Lt = Lt_mat.rows((j-1)*K*M,(j*K*M)-1);
        // inv_Ft = Ft_mat.rows((j-1)*M,(j*M)-1);

        rt_mat.col(j-1) = Z_t.t() * Ft_mat.rows((j-1)*M,(j*M)-1).t() * kal_resid.col(j-1) + arma::trans(Lt_mat.rows((j-1)*K*M,(j*K*M)-1))*rt_mat.col(j);
    }

    arma::mat beta_smoothed = arma::zeros(K*M,n_adj+1);

    for (int i=0; i < n_adj; i++) {
        beta_smoothed.col(i+1) = beta_smoothed.col(i) + Q_draw*rt_mat.col(i+1); // eq (8) in Durbin and Koopman
    }

    // beta_tilde = beta_smoothed + beta_plus;
    return beta_smoothed + beta_plus;
}
