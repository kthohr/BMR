/*################################################################################
  ##
  ##   Copyright (C) 2011-2017 Keith O'Hara
  ##
  ##   This file is part of the BMLib C++ library.
  ##
  ##   BMLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   BMLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/*
 * bvarw class
 */

#ifndef _bmlib_bvarw_HPP
#define _bmlib_bvarw_HPP

class bvarw
{
    public:
        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags
        int M;          // number of endogenous variables
        int K;          // number of coefficients in each 
        int n_ext_vars; // number of 'external' variables

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

        // member functions
        ~bvarw(){};
         bvarw(){};
        
        void build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp);
        void build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp);

        void prior(const arma::vec& coef_prior, const double Xi_beta, const arma::mat& Xi_Sigma, const int gamma);
        void prior(const arma::vec& coef_prior, const arma::mat& Xi_beta, const arma::mat& Xi_Sigma, const int gamma);

        void gibbs(const int n_draws, const int n_burnin);

        void IRF(const int n_irf_periods);

        arma::cube forecast(const int horizon, const bool incl_shocks);
        arma::cube forecast(const arma::mat& Y_T, const int horizon, const bool incl_shocks);

    protected:
        void build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp);
        arma::cube forecast_int(const arma::mat* Y_T_inp, const int horizon, const bool incl_shocks);
};

#endif
