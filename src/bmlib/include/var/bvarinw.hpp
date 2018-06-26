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
 * bvarinw class
 */

#ifndef _bmpp_bvarinw_HPP
#define _bmpp_bvarinw_HPP

class bvarinw
{
    public:
        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags
        int M;          // number of endogenous variables
        int K;          // number of coefficients in each
        int n_ext_vars; // number of 'external' variables

        bool only_stationary_draws = false; // discard non-stationary draws
        bool irfs_lr_restrict = false;      // impose long-run restictions on IRFs

        arma::mat Y; // Y = X beta + e
        arma::mat X;
        // arma::mat Z; // vec(Y) = Z alpha + vec(e)

        // ML-type estimates

        arma::mat alpha_hat;      // OLS estimate of alpha
        arma::mat Sigma_hat;      // OLS-based estimation of covariance matrix of 'e'

        // prior data

        arma::mat alpha_pr_mean;  // prior mean
        arma::mat alpha_pr_var;   // prior variance

        arma::mat Sigma_pr_scale; // prior scale matrix
        int Sigma_pr_dof;         // prior degrees of freedom

        // posterior data

        arma::mat alpha_pt_mean;  // posterior mean
        arma::mat alpha_pt_var;  // posterior mean

        int Sigma_pt_dof;         // posterior degrees of freedom
        arma::mat Sigma_pt_mean;  // posterior mean

        arma::cube beta_draws;    // posterior draws of beta
        arma::cube Sigma_draws;   // posterior draws of Sigma

        arma::cube irfs;          // irfs based on the posterior draws

        //
        // member functions

        ~bvarinw() = default;
         bvarinw() = default;

        bvarinw(const bvarinw&) = default;
        bvarinw& operator=(const bvarinw&) = default;

        bvarinw(bvarinw&&) = default;
        bvarinw& operator=(bvarinw&&) = default;

        void build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp);
        void build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp);

        void reset_draws();

        void prior(const arma::vec& coef_prior, const double Xi_beta, const double Xi_Sigma, const int gamma);
        void prior(const arma::vec& coef_prior, const arma::mat& Xi_beta, const arma::mat& Xi_Sigma, const int gamma);

        void gibbs(const uint_t n_draws, const uint_t n_burnin);

        arma::cube IRF(const uint_t n_irf_periods);
        arma::cube FEVD(const uint_t n_periods);

        arma::cube forecast(const uint_t horizon, const bool incl_shocks);
        arma::cube forecast(const arma::mat& Y_T, const uint_t horizon, const bool incl_shocks);

    protected:
        void build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp);
        arma::cube forecast_int(const arma::mat* Y_T_inp, const uint_t horizon, const bool incl_shocks);
};

#endif
