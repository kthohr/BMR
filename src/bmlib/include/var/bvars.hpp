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
 * bvars class
 */

#ifndef _bmpp_bvars_HPP
#define _bmpp_bvars_HPP

class bvars
{
    public:
        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags; 'k' in Villani's notation
        int q;          // number of external variables + c_int
        int M;          // number of endogenous variables; 'p' in Villani's notation
        int K;          // M*p
        int n_ext_vars; // number of 'external' variables

        bool only_stationary_draws = false; // discard non-stationary draws
        bool irfs_lr_restrict = false;      // impose long-run restictions on IRFs

        arma::mat Y;   // beta(L) (Y_t - Psi*d_t) = e_t
        arma::mat X;   // Lags of Y

        arma::mat d;   // [d_t]
        arma::mat d_X; // [d_{t-1},...,d_{t-p}]

        arma::mat D;   // 'exogeneous' variables (including a constant), D = [d_t, - d_{t-1}, ..., -d_{t-p}]

        // ML-type estimates

        arma::mat psi_hat;
        arma::mat alpha_hat;       // OLS estimate of beta
        arma::mat Sigma_hat;      // OLS-based estimation of covariance matrix of 'e'

        // prior data

        arma::mat psi_pr_mean;    // psi = vec(Psi)
        arma::mat psi_pr_var;

        arma::mat alpha_pr_mean;  // prior mean of alpha
        arma::mat alpha_pr_var;   // prior variance of alpha

        arma::mat Sigma_pr_scale; // prior scale matrix
        int Sigma_pr_dof;         // prior degrees of freedom

        // posterior data

        arma::vec psi_pt_mean;    // posterior mean of psi = vec(Psi)
        arma::mat psi_pt_var;     // posterior variance of psi

        arma::vec alpha_pt_mean;  // posterior mean of alpha
        arma::mat alpha_pt_var;   // posterior variance of alpha

        int Sigma_pt_dof;         // posterior degrees of freedom
        arma::mat Sigma_pt_mean;  // posterior mean

        arma::cube Psi_draws;     // posterior draws of Psi
        arma::cube beta_draws;    // posterior draws of beta
        arma::cube Sigma_draws;   // posterior draws of Sigma

        //
        // member functions

        ~bvars() = default;
         bvars() = default;

        bvars(const bvars&) = default;
        bvars& operator=(const bvars&) = default;

        bvars(bvars&&) = default;
        bvars& operator=(bvars&&) = default;
        
        void build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp);
        void build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp);

        void reset_draws();

        void prior(const arma::vec& coef_prior, const double HP_1, const double HP_2, 
                   const arma::mat& Psi_prior, const double Xi_psi, const int gamma, 
                   const bool full_cov_prior = true);

        void gibbs(const uint_t n_draws, const uint_t n_burnin);

        arma::cube IRF(const uint_t n_irf_periods);
        arma::cube FEVD(const uint_t n_periods);

        arma::cube forecast(const uint_t horizon, const bool incl_shocks);
        arma::cube forecast(const arma::mat& Y_T, const uint_t horizon, const bool incl_shocks);

    private:
        void build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp);
        arma::mat minn_pr_var();
        arma::cube forecast_int(const arma::mat* Y_T_inp, const uint_t horizon, const bool incl_shocks);
};

#endif
