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
 * DSGE-VAR model class
 */

#ifndef _bmpp_dsgevar_class_HPP
#define _bmpp_dsgevar_class_HPP

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

        double lambda = 1.0; // weight on DSGE prior

        // a DSGE object

        mutable dsge<T> dsge_obj;

        // data used for estimation

        arma::mat Y;
        arma::mat X;

        arma::mat YY;
        arma::mat XY;
        arma::mat XX;

        // posterior data

        arma::vec alpha_pt_mean;   // posterior mean
        arma::mat alpha_pt_var;    // posterior variance

        int Sigma_pt_dof;         // posterior degrees of freedom
        arma::mat Sigma_pt_mean;  // posterior mean

        arma::cube beta_draws;    // posterior draws of beta
        arma::cube Sigma_draws;   // posterior draws of Sigma

        //
        // member functions

        ~dsgevar() = default;
         dsgevar() = default;

        dsgevar(const dsgevar&) = default;
        dsgevar& operator=(const dsgevar&) = default;

        dsgevar(dsgevar&&) = default;
        dsgevar& operator=(dsgevar&&) = default;

        void build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp, const double lambda_inp);

        void set_bounds(const arma::vec& lower_bounds_inp, const arma::vec& upper_bounds_inp);
        void set_prior(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp);

        double log_posterior_kernel(const arma::vec& pars_inp) const;

        arma::vec estim_mode(const arma::vec& initial_vals);
        arma::vec estim_mode(const arma::vec& initial_vals, arma::mat& vcov_mat);
        arma::vec estim_mode(const arma::vec& initial_vals, arma::mat* vcov_mat, optim::algo_settings* settings_inp);

        void estim_mcmc(const arma::vec& initial_vals, mcmc::algo_settings_t* settings_inp);

        arma::cube IRF(const int n_irf_periods) const;

        arma::cube forecast(const int n_horizon, const bool incl_shocks);
        arma::cube forecast(const arma::mat& X_T, const int n_horizon, const bool incl_shocks);

        arma::cube state_filter() const;

    protected:
        void model_moments(arma::mat& Gamma_YY, arma::mat& Gamma_XY, arma::mat& Gamma_XX) const;
        void model_moments(const dsge<T>& dsge_model_inp, arma::mat& Gamma_YY, arma::mat& Gamma_XY, arma::mat& Gamma_XX) const;

        double log_likelihood(const arma::mat& Gamma_YY, const arma::mat& Gamma_XY, const arma::mat& Gamma_XX, const arma::mat& Gamma_bar_YY, const arma::mat& Gamma_bar_XY, const arma::mat& Gamma_bar_XX) const;
        double log_likelihood_inf(const arma::mat& Gamma_YY, const arma::mat& Gamma_XY, const arma::mat& Gamma_XX) const;

        void gibbs();

        static double mode_objfn(const arma::vec& pars_inp, arma::vec* grad_vec, void* mode_data);
        static double mcmc_objfn(const arma::vec& pars_inp, void* mode_data);

        arma::cube forecast_int(const arma::mat* X_T_inp, const int horizon, const bool incl_shocks);
};

template<typename T>
struct dsgevar_estim_data {
    dsgevar<T> dsgevar_obj;
};

#include "dsgevar_class.tpp"

#endif
