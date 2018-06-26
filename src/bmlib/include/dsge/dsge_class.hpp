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
 * DSGE model class
 */

#ifndef _bmpp_dsge_class_HPP
#define _bmpp_dsge_class_HPP

template<class T>
class dsge
{
    public:

        // a LREM object (e.g., gensys or uhlig)

        mutable T lrem_obj;

        // state equation:       X_t = F*X_{t-1} + G*shock_t,
        // measurement equation: Y_t = C + H'*X_t + measerr_t

        mutable arma::mat kalman_mat_F; // state VAR(1) matrix
        mutable arma::mat kalman_mat_Q; // covariance matrix = SS_G * cov_shock_t * SS_G.t()

        mutable arma::mat kalman_mat_C;
        mutable arma::mat kalman_mat_H;
        mutable arma::mat kalman_mat_R; // covariance matrix of meas. err.

        // priors and bounds

        arma::uvec prior_form;
        arma::mat prior_pars;

        arma::vec lower_bounds;
        arma::vec upper_bounds;

        // choice of filtering method

        int filter_choice = 1; // 1 for Kalman, 2 for Chandrasekhar

        // data used for estimation

        arma::mat estim_data;

        // MCMC results

        arma::mat dsge_draws;

        // a function mapping the 'deep' parameters to the structural matrices

        std::function<void (const arma::vec& pars_inp, T& lrem_obj_inp, arma::mat& shocks_cov_out, arma::mat& C_out, arma::mat& H_out, arma::mat& R_out)> model_fn;

        //
        // member functions

        ~dsge() = default;
         dsge() = default;

        dsge(const dsge&) = default;
        dsge& operator=(const dsge&) = default;

        dsge(dsge&&) = default;
        dsge& operator=(dsge&&) = default;

        void set_bounds(const arma::vec& lower_bounds_inp, const arma::vec& upper_bounds_inp);
        void set_prior(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp);

        void solve() const;
        void state_space(arma::mat& F_state, arma::mat& G_state) const;
        arma::mat simulate(const int sim_periods, const int burnin) const;

        void solve_to_state_space(const arma::vec& pars_inp) const;

        double log_prior(const arma::vec& pars_inp) const;
        double log_posterior_kernel(const arma::vec& pars_inp) const;

        arma::vec estim_mode(const arma::vec& initial_vals);
        arma::vec estim_mode(const arma::vec& initial_vals, arma::mat& vcov_mat);
        arma::vec estim_mode(const arma::vec& initial_vals, arma::mat* vcov_mat, optim::algo_settings* settings_inp);

        void estim_mcmc(const arma::vec& initial_vals);
        void estim_mcmc(const arma::vec& initial_vals, mcmc::algo_settings_t* settings_inp);

        arma::cube IRF(const int n_irf_periods, const bool observ_irfs) const;

        arma::cube forecast(const int n_horizon, const bool incl_shocks) const;

        arma::cube state_filter() const;

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
