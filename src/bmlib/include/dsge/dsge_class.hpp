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
 * DSGE model class
 */

#ifndef _bmlib_dsge_HPP
#define _bmlib_dsge_HPP

template<class T>
class dsge
{
    public:

        // a LREM object (e.g., gensys or uhlig)

        T lrem_obj;

        // state equation:       X_t = F*X_{t-1} + G*shock_t,
        // measurement equation: Y_t = C + H'*X_t + measerr_t

        arma::mat kalman_mat_F; // state VAR(1) matrix
        arma::mat kalman_mat_Q; // covariance matrix = SS_G * cov_shock_t * SS_G.t()

        arma::mat kalman_mat_C;
        arma::mat kalman_mat_H;
        arma::mat kalman_mat_R; // covariance matrix of meas. err.

        // priors and bounds

        arma::uvec prior_form;
        arma::mat prior_pars;

        arma::vec lower_bounds;
        arma::vec upper_bounds;

        // data used for estimation

        arma::mat estim_data;

        // MCMC results

        arma::mat mcmc_dsge_pars;

        // a function mapping the 'deep' parameters to the structural matrices

        std::function<void (const arma::vec& pars_inp, T& lrem_obj, arma::mat& Q_out, arma::mat& C_out, arma::mat& R_out)> pars_to_mats;

        //
        // member functions

        void set_bounds(const arma::vec& lower_bounds_inp, const arma::vec& upper_bounds_inp);
        void set_prior(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp);

        void solve();
        void state_space(arma::mat& F_state, arma::mat& G_state);
        arma::mat simulate(const int sim_periods, const int burnin);

        void solve_to_state_space(const arma::vec& pars_inp);

        double log_prior(const arma::vec& pars_inp);
        double log_posterior_kernel(const arma::vec& pars_inp);

        arma::vec estim_mode(const arma::vec& initial_vals);
        arma::vec estim_mode(const arma::vec& initial_vals, optim::opt_settings* settings_inp);

        void estim_mcmc(const arma::vec& initial_vals);
        void estim_mcmc(const arma::vec& initial_vals, mcmc::mcmc_settings* settings_inp);

    protected:
        static double mode_objfn(const arma::vec& pars_inp, arma::vec* grad_vec, void* mode_data);
        static double mcmc_objfn(const arma::vec& pars_inp, void* mode_data);
};

template<typename T>
struct dsge_estim_data {
    dsge<T> dsge_obj;
};

#include "dsge_class.tpp"

#endif
