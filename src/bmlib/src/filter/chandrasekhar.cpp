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
 * Chandrasekhar recursions
 */

#include "misc/misc.hpp"
namespace bm {
    #include "filter/chandrasekhar.hpp"
}

// X_t = F*X_{t-1} + e, e ~ N(0,Q)
// Y_t = H'*X_{t}  + z, z ~ N(0,R)

double 
bm::chand_recur(const arma::mat& data_inp, const arma::mat& F, const arma::mat& Q, const arma::mat& C, const arma::mat& H, const arma::mat& R)
{
    //
    const int n = data_inp.n_rows;
    const int k = data_inp.n_cols;
    
    const double nlog2pi = k*log(2*arma::datum::pi);
    
    const arma::mat tH = H.t();
    const arma::mat tF = F.t();

    // Solve for the steady-state covariance matrix of the state vector using a doubling algorithm

    arma::mat state_pos_filt = arma::zeros(F.n_cols,1);
    
    arma::mat state_cov_filt = lyapunov_dbl(Q, F);
    
    // Now initialize the filter using the steady-state covariance matrix

    // arma::mat state_pos_filt = InitialState;
    arma::mat state_cov_pred = F*state_cov_filt*tF + Q;

    arma::mat Sigma = tH*state_cov_pred*H + R;
    arma::mat inv_Sigma = arma::inv_sympd(Sigma);
    
    arma::mat St = F*state_cov_pred*H;
    arma::mat Mt = -inv_Sigma;    
    arma::mat Kt = St*inv_Sigma;
    
    arma::mat resid = data_inp.row(0).t() - tH*state_pos_filt - C;
    
    double loglikelihood = nlog2pi + std::log(arma::det(Sigma)) + arma::as_scalar(resid.t()*inv_Sigma*resid);
    
    //
    // begin loop

    arma::mat Sigma_p = Sigma, inv_Sigma_p = inv_Sigma;

    for (int i=1; i < n; i++) {
        state_pos_filt = F*state_pos_filt + Kt*resid;

        resid = data_inp.row(i).t() - tH*state_pos_filt - C;
        
        arma::mat tHSt = tH*St;
        arma::mat MSpZp = Mt*tHSt.t();
        arma::mat FSt = F*St;
        
        Sigma_p = Sigma + tHSt*MSpZp;     
        inv_Sigma_p = arma::inv_sympd(Sigma_p);

        //

        Kt = (Kt*Sigma + FSt*MSpZp)*inv_Sigma_p;
        St = FSt - Kt*tHSt;
        Mt += MSpZp * inv_Sigma * MSpZp.t();
    
        Sigma = Sigma_p;
        inv_Sigma = inv_Sigma_p;
        
        loglikelihood += nlog2pi + std::log(arma::det(Sigma)) + arma::as_scalar(resid.t()*inv_Sigma*resid);
    }
    
    // Returns the loglikelihood value
    return -0.5*loglikelihood;
}
