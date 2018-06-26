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
 * Kalman filter and smoother
 */

#include "misc/misc.hpp"
namespace bm {
    #include "filter/kalman.hpp"
}

// X_t = F*X_{t-1} + e, e ~ N(0,Q)
// Y_t = H'*X_{t}  + z, z ~ N(0,R)

double
bm::kalman_filter(const arma::mat& data_inp, const arma::mat& F, const arma::mat& Q, const arma::mat& C, const arma::mat& H, const arma::mat& R, arma::mat* state_filt_out)
{
    //
    const int n = data_inp.n_rows;
    const int k = data_inp.n_cols;
    
    const double nlog2pi = k*log(2*arma::datum::pi);
    
    const arma::mat tH = H.t();
    const arma::mat tF = F.t();

    if (state_filt_out) {
        state_filt_out->set_size(n, (int) F.n_cols);
    }

    // Solve for the steady-state covariance matrix of the state vector using a doubling algorithm

    arma::mat state_pos_filt = arma::zeros(F.n_cols,1);
    
    arma::mat state_cov_filt = lyapunov_dbl(Q, F);
    
    // Now initialize the filter using the steady-state covariance matrix

    arma::mat state_pos_pred = F*state_pos_filt;
    arma::mat state_cov_pred = F*state_cov_filt*tF + Q;

    arma::mat Sigma = tH*state_cov_pred*H + R;
    arma::mat inv_Sigma = arma::inv_sympd(Sigma);
    //
    arma::mat kal_gain = state_cov_pred*H*inv_Sigma;
    arma::mat resid = data_inp.row(0).t() - tH*state_pos_pred - C;

    state_pos_filt = state_pos_pred + kal_gain*resid;
    state_cov_filt = state_cov_pred - kal_gain*tH*state_cov_pred;

    double loglikelihood = nlog2pi + std::log(arma::det(Sigma)) + arma::as_scalar(resid.t()*inv_Sigma*resid);

    if (state_filt_out) {
        state_filt_out->row(0) = state_pos_filt.t();
    }
    
    //
    // begin loop

    for (int i=1; i < n; i++) {
        state_pos_pred = F*state_pos_filt;
        state_cov_pred = F*state_cov_filt*tF + Q;

        Sigma = tH*state_cov_pred*H + R;
        inv_Sigma = arma::inv_sympd(Sigma);

        //

        kal_gain = state_cov_pred*H*inv_Sigma;
        resid = data_inp.row(i).t() - tH*state_pos_pred - C;

        state_pos_filt = state_pos_pred + kal_gain*resid;
        state_cov_filt = state_cov_pred - kal_gain*tH*state_cov_pred;

        loglikelihood += nlog2pi + std::log(arma::det(Sigma)) + arma::as_scalar(resid.t()*inv_Sigma*resid);

        if (state_filt_out) {
            state_filt_out->row(i) = state_pos_filt.t();
        }
    }

    // Returns the loglikelihood value
    return -0.5*loglikelihood;
}

double
bm::kalman_filter(const arma::mat& data_inp, const arma::mat& F, const arma::mat& Q, const arma::mat& C, const arma::mat& H, const arma::mat& R)
{
    return kalman_filter(data_inp,F,Q,C,H,R,nullptr);
}
