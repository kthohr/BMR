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
 * bvartvp class
 */

#ifndef _bmlib_bvartvp_HPP
#define _bmlib_bvartvp_HPP

class bvartvp
{
    public:
        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags
        int M;          // number of endogenous variables
        int K;          // number of coefficients in each
        int n_ext_vars; // number of 'external' variables

        int tau;        // training sample length

        arma::mat Y;
        arma::mat Z;

        // ML-type estimates
        arma::mat alpha_hat;      // OLS estimate of alpha
        arma::mat Sigma_hat;      // OLS-based estimation of covariance matrix of 'e'

        // prior data
        arma::mat alpha_pr_mean;  // prior mean
        arma::mat alpha_pr_var;   // prior variance

        arma::mat Q_pr_scale;     // prior scale matrix
        int Q_pr_dof;             // prior degrees of freedom

        arma::mat Sigma_pr_scale; // prior scale matrix
        int Sigma_pr_dof;         // prior degrees of freedom

        // posterior data
        arma::mat alpha_pt_mean;  // posterior mean

        int Q_pt_dof;             // posterior degrees of freedom
        arma::mat Q_pt_mean;      // posterior mean

        int Sigma_pt_dof;         // posterior degrees of freedom
        arma::mat Sigma_pt_mean;  // posterior mean

        arma::cube alpha_draws;   // posterior draws of alpha
        arma::cube Q_draws;       // posterior draws of Q
        arma::cube Sigma_draws;   // posterior draws of Sigma

        arma::cube irfs;          // irfs based on the posterior draws

        // member functions
        ~bvartvp(){};
         bvartvp(){};

        void build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp);
        // void build(const arma::mat& data_raw, const arma::mat& data_ext);

        void reset_draws();

        void prior(const int tau_inp, const double Xi_beta, const double Xi_Q, const int gamma_Q, const double Xi_Sigma, const int gamma_S);

        void gibbs(const int n_draws, const int n_burnin);

        void IRF(const int n_irf_periods, const int time_ind);

        arma::cube forecast(const int horizon, const bool incl_shocks);
        arma::cube forecast(const arma::mat& Y_T, const int horizon, const bool incl_shocks);

    private:
        void build_int(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp);
        arma::cube forecast_int(const arma::mat* Y_T_inp, const int horizon, const bool incl_shocks);
};

#endif
