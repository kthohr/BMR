/*################################################################################
  ##
  ##   Copyright (C) 2016-2017 Keith O'Hara
  ##
  ##   This file is part of the OptimLib C++ library.
  ##
  ##   OptimLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   OptimLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/*
 * Particle Swarm Optimization (PSO) with Differentially-Perturbed Velocity (DV)
 *
 * Keith O'Hara
 * 08/04/2016
 *
 * This version:
 * 08/14/2017
 */

#include "optim.hpp"

bool
optim::pso_dv_int(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)> opt_objfn, void* opt_data, opt_settings* settings_inp)
{
    bool success = false;

    const double BIG_POS_VAL = OPTIM_BIG_POS_VAL;
    const int n_vals = init_out_vals.n_elem;

    //
    //

    opt_settings settings;

    if (settings_inp) {
        settings = *settings_inp;
    }

    const int conv_failure_switch = settings.conv_failure_switch;
    const double err_tol = settings.err_tol;

    const int n_pop = (settings.pso_n_pop > 0) ? settings.pso_n_pop : 100;
    const int n_gen = (settings.pso_n_gen > 0) ? settings.pso_n_gen : 1000;

    const int stag_limit = 50;

    double par_w = 1.0;
    double par_beta = 0.5;
    const double par_damp = 0.99;
    // const double par_c_1 = 1.494;
    const double par_c_2 = 1.494;

    const double par_CR = 0.7;

    const arma::vec par_initial_lb = ((int) settings.pso_initial_lb.n_elem == n_vals) ? settings.pso_initial_lb : init_out_vals - 0.5;
    const arma::vec par_initial_ub = ((int) settings.pso_initial_ub.n_elem == n_vals) ? settings.pso_initial_ub : init_out_vals + 0.5;

    const bool vals_bound = settings.vals_bound;
    
    const arma::vec lower_bounds = settings.lower_bounds;
    const arma::vec upper_bounds = settings.upper_bounds;

    const arma::uvec bounds_type = determine_bounds_type(vals_bound, n_vals, lower_bounds, upper_bounds);

    // lambda function for box constraints

    std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* box_data)> box_objfn = [opt_objfn, vals_bound, bounds_type, lower_bounds, upper_bounds] (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data) -> double {
        //

        if (vals_bound) {
            arma::vec vals_inv_trans = inv_transform(vals_inp, bounds_type, lower_bounds, upper_bounds);
            
            return opt_objfn(vals_inv_trans,nullptr,opt_data);
        } else {
            return opt_objfn(vals_inp,nullptr,opt_data);
        }
    };

    //
    // initialize

    arma::vec objfn_vals(n_pop);
    arma::mat P(n_pop,n_vals);

#ifdef OPTIM_OMP
    #pragma omp parallel for
#endif
    for (int i=0; i < n_pop; i++) {
        P.row(i) = par_initial_lb.t() + (par_initial_ub.t() - par_initial_lb.t())%arma::randu(1,n_vals);

        double prop_objfn_val = opt_objfn(P.row(i).t(),nullptr,opt_data);

        if (std::isnan(prop_objfn_val)) {
            prop_objfn_val = BIG_POS_VAL;
        }
        
        objfn_vals(i) = prop_objfn_val;

        if (vals_bound) {
            P.row(i) = arma::trans( transform(P.row(i).t(), bounds_type, lower_bounds, upper_bounds) );
        }
    }

    arma::vec best_vals = objfn_vals;

    arma::mat best_vecs = P;

    arma::mat V = arma::zeros(n_pop,n_vals);

    double global_best_val = objfn_vals.min();
    arma::rowvec global_best_vec = P.row( objfn_vals.index_min() );

    arma::vec stag_vec = arma::zeros(n_pop,1);

    //
    //

    int iter = 0;
    double err = 2.0*err_tol;

    while (err > err_tol && iter < n_gen) {
        iter++;

        arma::rowvec P_max = arma::max(P);
        arma::rowvec P_min = arma::min(P);

#ifdef OPTIM_OMP
        #pragma omp parallel for 
#endif
        for (int i=0; i < n_pop; i++) {

            int c_1, c_2;

            do { // 'r_2' in paper's notation
                c_1 = arma::as_scalar(arma::randi(1, arma::distr_param(0, n_pop-1)));
            } while(c_1==i);

            do { // 'r_3' in paper's notation
                c_2 = arma::as_scalar(arma::randi(1, arma::distr_param(0, n_pop-1)));
            } while(c_2==i || c_2==c_1);

            //

            arma::vec rand_CR = arma::randu(n_vals,1);

            arma::rowvec delta_vec = P.row(c_1) - P.row(c_2);

            for (int k=0; k < n_vals; k++) {
                if (rand_CR(k) <= par_CR) {
                    double rand_u = arma::as_scalar(arma::randu(1));
                    V(i,k) = par_w*V(i,k) + par_beta*delta_vec(k) + par_c_2*rand_u*(global_best_vec(k) - P(i,k));
                }
            }

            arma::rowvec TR = P.row(i) + V.row(i);
            double TR_objfn_val = box_objfn(TR.t(),nullptr,opt_data);

            if (TR_objfn_val < objfn_vals(i)) {
                P.row(i) = TR;
                objfn_vals(i) = TR_objfn_val;
            } else {
                stag_vec(i) += 1;
            }

            if (stag_vec(i) >= stag_limit) {
                P.row(i) = P_min + arma::randu(1,n_vals) % (P_max - P_min);
                stag_vec(i) = 0;

                objfn_vals(i) = box_objfn(P.row(i).t(),nullptr,opt_data);
            }
                
            // if (objfn_vals(i) < best_vals(i)) {
            //     best_vals(i) = objfn_vals(i);
            //     best_vecs.row(i) = P.row(i);
            // }
        }

        if (objfn_vals.min() < global_best_val) {
            global_best_val = objfn_vals.min();
            global_best_vec = P.row( objfn_vals.index_min() );
        }

        par_w *= par_damp;
        // par_w = std::min(0.4,par_w*par_damp);
    }
    //
    if (vals_bound) {
	    global_best_vec = arma::trans( inv_transform(global_best_vec.t(), bounds_type, lower_bounds, upper_bounds) );
    }

    error_reporting(init_out_vals,global_best_vec.t(),opt_objfn,opt_data,success,err,err_tol,iter,n_gen,conv_failure_switch,settings_inp);
    //
    return true;
}

bool
optim::pso_dv(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)> opt_objfn, void* opt_data)
{
    return pso_dv_int(init_out_vals,opt_objfn,opt_data,nullptr);
}

bool
optim::pso_dv(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)> opt_objfn, void* opt_data, opt_settings& settings)
{
    return pso_dv_int(init_out_vals,opt_objfn,opt_data,&settings);
}
