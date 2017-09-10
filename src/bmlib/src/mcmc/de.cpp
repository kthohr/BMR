/*################################################################################
  ##
  ##   Copyright (C) 2011-2017 Keith O'Hara
  ##
  ##   This file is part of the MCMC C++ library.
  ##
  ##   MCMC is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   MCMC is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/
 
/*
 * Differential Evolution (DE) MCMC
 *
 * Keith O'Hara
 * 03/01/2016
 *
 * This version:
 * 08/12/2017
 */

#include "mcmc.hpp"

bool
mcmc::de_int(const arma::vec& initial_vals, arma::cube& draws_out, std::function<double (const arma::vec& vals_inp, void* target_data)> target_log_kernel, void* target_data, mcmc_settings* settings_inp)
{
    bool success = false;

    const double BIG_NEG_VAL = MCMC_BIG_NEG_VAL;
	const int n_vals = initial_vals.n_elem;
    
    //
    // DE settings

    mcmc_settings settings;

    if (settings_inp) {
        settings = *settings_inp;
    }

    const int n_pop    = settings.de_n_pop;
    const int n_gen    = settings.de_n_gen;
    const int n_burnin = settings.de_n_burnin;

    const bool jumps = settings.de_jumps;
    const double par_b = settings.de_par_b;
    // const double par_gamma = settings.de_par_gamma;
    const double par_gamma = 2.38 / std::sqrt(2.0*n_vals);
    const double par_gamma_jump = settings.de_par_gamma_jump;

    const arma::vec par_initial_lb = ((int) settings.de_initial_lb.n_elem == n_vals) ? settings.de_initial_lb : initial_vals - 0.5;
    const arma::vec par_initial_ub = ((int) settings.de_initial_ub.n_elem == n_vals) ? settings.de_initial_ub : initial_vals + 0.5;

    const bool vals_bound = settings.vals_bound;

    const arma::vec lower_bounds = settings.lower_bounds;
    const arma::vec upper_bounds = settings.upper_bounds;

    const arma::uvec bounds_type = determine_bounds_type(vals_bound, n_vals, lower_bounds, upper_bounds);

    // lambda function for box constraints

    std::function<double (const arma::vec& vals_inp, void* box_data)> box_log_kernel = [target_log_kernel, vals_bound, bounds_type, lower_bounds, upper_bounds] (const arma::vec& vals_inp, void* target_data) -> double {
        //
        if (vals_bound) {
            arma::vec vals_inv_trans = inv_transform(vals_inp, bounds_type, lower_bounds, upper_bounds);

            return target_log_kernel(vals_inv_trans, target_data) + log_jacobian(vals_inp, bounds_type, lower_bounds, upper_bounds);
        } else {
            return target_log_kernel(vals_inp, target_data);
        }
    };

    //
    arma::vec target_vals(n_pop);
    arma::mat X(n_pop,n_vals);

#ifdef MCMC_OMP
    #pragma omp parallel for
#endif
    for (int i=0; i < n_pop; i++) {
        X.row(i) = par_initial_lb.t() + (par_initial_ub.t() - par_initial_lb.t())%arma::randu(1,n_vals);

        double prop_kernel_val = box_log_kernel(X.row(i).t(),target_data);

        if (!std::isfinite(prop_kernel_val)) {
            prop_kernel_val = BIG_NEG_VAL;
        }
        
        target_vals(i) = prop_kernel_val;
    }

    //
    // begin MCMC run

    draws_out.set_size(n_pop,n_vals,n_gen);

    int n_accept = 0;
    double par_gamma_run = par_gamma;
    
    for (int j=0; j < n_gen + n_burnin; j++) {
        double temperature_j = de_cooling_schedule(j,n_gen);
        
        if (jumps && ((j+1) % 10 == 0)) {
            par_gamma_run = par_gamma_jump;
        }

#ifdef MCMC_OMP
        #pragma omp parallel for
#endif
            for (int i=0; i < n_pop; i++) {

                int R_1, R_2;

                do {
                    R_1 = arma::as_scalar(arma::randi(1, arma::distr_param(0, n_pop-1)));
                } while(R_1==i);

                do {
                    R_2 = arma::as_scalar(arma::randi(1, arma::distr_param(0, n_pop-1)));
                } while(R_2==i || R_2==R_1);

                //

                arma::vec prop_rand = arma::randu(n_vals,1)*(2*par_b) - par_b; // generate a vector of U[-b,b] RVs

                arma::rowvec X_prop = X.row(i) + par_gamma_run * ( X.row(R_1) - X.row(R_2) ) + prop_rand;
                
                double prop_kernel_val = box_log_kernel(X_prop.t(),target_data);
                
                if (!std::isfinite(prop_kernel_val)) {
                    prop_kernel_val = BIG_NEG_VAL;
                }
                
                double comp_val = prop_kernel_val - target_vals(i);
                //
                if (comp_val > temperature_j * std::log(arma::as_scalar(arma::randu(1)))) {
                    X.row(i) = X_prop;
                    
                    target_vals(i) = prop_kernel_val;
                    
                    if (j >= n_burnin) {
                        n_accept++;
                    }
                }
            }
        //
        if(j >= n_burnin){
            draws_out.slice(j-n_burnin) = X;
        }
        //
        if (jumps && ((j+1) % 10 == 0)) {
            par_gamma_run = par_gamma;
        }
    }
	//
	if (vals_bound) {
#ifdef MCMC_OMP
        #pragma omp parallel for
#endif
        for (int ii = 0; ii < n_gen; ii++) {
            for (int jj = 0; jj < n_pop; jj++) {
                draws_out.slice(ii).row(jj) = arma::trans(inv_transform(draws_out.slice(ii).row(jj).t(), bounds_type, lower_bounds, upper_bounds));
            }
        }
	}
	//
    if (settings_inp) {
	    settings_inp->de_accept_rate = (double) n_accept / (double) (n_pop*n_gen);
    }
	//
    success = true;
	return success;
}

// wrappers

bool
mcmc::de(const arma::vec& initial_vals, arma::cube& draws_out, std::function<double (const arma::vec& vals_inp, void* target_data)> target_log_kernel, void* target_data)
{
    return de_int(initial_vals,draws_out,target_log_kernel,target_data,nullptr);
}

bool
mcmc::de(const arma::vec& initial_vals, arma::cube& draws_out, std::function<double (const arma::vec& vals_inp, void* target_data)> target_log_kernel, void* target_data, mcmc_settings& settings)
{
    return de_int(initial_vals,draws_out,target_log_kernel,target_data,&settings);
}
