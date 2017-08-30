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
 * DSGE-VAR model class
 */

#ifndef _bmlib_dsgevar_HPP
#define _bmlib_dsgevar_HPP

template<class T>
class dsgevar
{
    public:
        //
        // objects

        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags
        int K;          // number of coefficients in each
        int M;          // number of endogenous variables

        bool finite_lambda;
        double lambda = 1.0; // weight on DSGE prior

        // a DSGE object

        dsge<T> dsge_obj;

        // data used for estimation

        arma::mat Y;
        arma::mat X;

        arma::mat YY;
        arma::mat XY;
        arma::mat XX;

        // posterior data

        arma::vec beta_pt_mean;  // posterior mean
        arma::mat beta_pt_var;   // posterior variance

        int Sigma_pt_dof;         // posterior degrees of freedom
        arma::mat Sigma_pt_mean;  // posterior mean

        arma::cube beta_draws;    // posterior draws of beta
        arma::cube Sigma_draws;   // posterior draws of Sigma
        
        arma::cube irfs;          // irfs based on the posterior draws

        //
        // member functions

        void build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp, const double lambda_inp);

        void set_bounds(const arma::vec& lower_bounds_inp, const arma::vec& upper_bounds_inp);
        void set_prior(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp);

        void model_moments(arma::mat& Gamma_YY, arma::mat& Gamma_XY, arma::mat& Gamma_XX);
        void model_moments(const dsge<T>& dsge_model_inp, arma::mat& Gamma_YY, arma::mat& Gamma_XY, arma::mat& Gamma_XX);

        double log_likelihood(const arma::mat& Gamma_YY, const arma::mat& Gamma_XY, const arma::mat& Gamma_XX, const arma::mat& Gamma_bar_YY, const arma::mat& Gamma_bar_XY, const arma::mat& Gamma_bar_XX);
        double log_likelihood_inf(const arma::mat& Gamma_YY, const arma::mat& Gamma_XY, const arma::mat& Gamma_XX);

        double log_posterior_kernel(const arma::vec& pars_inp);

        arma::vec estim_mode(const arma::vec& initial_vals);
        arma::vec estim_mode(const arma::vec& initial_vals, optim::opt_settings* settings_inp);

        void gibbs();
        void estim_mcmc(const arma::vec& initial_vals, mcmc::mcmc_settings* settings_inp);

        void IRF(const int n_irf_periods);

    protected:
        static double mode_objfn(const arma::vec& pars_inp, arma::vec* grad_vec, void* mode_data);
        static double mcmc_objfn(const arma::vec& pars_inp, void* mode_data);
};

template<typename T>
struct dsgevar_estim_data {
    dsgevar<T> dsgevar_obj;
};

#include "dsgevar_class.tpp"

#endif
