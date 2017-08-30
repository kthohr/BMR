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
 * BFGS method for quasi-Newton-based non-linear optimization
 *
 * Keith O'Hara
 * 12/23/2016
 *
 * This version:
 * 08/14/2017
 */

#include "optim.hpp"

bool
optim::bfgs_int(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)> opt_objfn, void* opt_data, opt_settings* settings_inp)
{   // notation: 'p' stands for '+1'.
    //
    bool success = false;

    const int n_vals = init_out_vals.n_elem;

    //
    // BFGS settings

    opt_settings settings;

    if (settings_inp) {
        settings = *settings_inp;
    }
    
    const int conv_failure_switch = settings.conv_failure_switch;
    const int iter_max = settings.iter_max;
    const double err_tol = settings.err_tol;

    const double wolfe_cons_1 = 1E-03; // line search tuning parameters
    const double wolfe_cons_2 = 0.90;

    const bool vals_bound = settings.vals_bound;
    
    const arma::vec lower_bounds = settings.lower_bounds;
    const arma::vec upper_bounds = settings.upper_bounds;

    const arma::uvec bounds_type = determine_bounds_type(vals_bound, n_vals, lower_bounds, upper_bounds);

    // lambda function for box constraints

    std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* box_data)> box_objfn = [opt_objfn, vals_bound, bounds_type, lower_bounds, upper_bounds] (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data) -> double {
        //

        if (vals_bound) {

            arma::vec vals_inv_trans = inv_transform(vals_inp, bounds_type, lower_bounds, upper_bounds);
            
            double ret;
            
            if (grad_out) {
                arma::vec grad_obj = *grad_out;

                ret = opt_objfn(vals_inv_trans,&grad_obj,opt_data);

                arma::mat jacob_matrix = jacobian_adjust(vals_inp,bounds_type,lower_bounds,upper_bounds);

                // *grad_out = jacob_matrix.t() * grad_obj; // correct gradient for transformation
                *grad_out = jacob_matrix * grad_obj; // no need for transpose as jacob_matrix is diagonal
            } else {
                ret = opt_objfn(vals_inv_trans,nullptr,opt_data);
            }

            return ret;
        } else {
            double ret = opt_objfn(vals_inp,grad_out,opt_data);

            return ret;
        }
    };

    //
    // initialization

    arma::vec x = init_out_vals;

    if (!x.is_finite()) {
        printf("bfgs error: non-finite initial value(s).\n");
        
        return false;
    }

    if (vals_bound) { // should we transform the parameters?
	    x = transform(x, bounds_type, lower_bounds, upper_bounds);
    }

    arma::mat W = arma::eye(n_vals,n_vals); // initial approx. to (inverse) Hessian 
    const arma::mat I_mat = arma::eye(n_vals,n_vals);

    arma::vec grad(n_vals); // gradient vector
    box_objfn(x,&grad,opt_data);

    // double err = arma::accu(arma::abs(grad));
    double err = arma::norm(grad, 2);
    if (err <= err_tol) {
        return true;
    }

    //
    // if ||gradient(initial values)|| > tolerance, then continue

    // double t_line = 1.0;    // initial line search value
    arma::vec d = - W*grad; // direction

    arma::vec x_p = x, grad_p = grad;

    line_search_mt(1.0, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, box_objfn, opt_data);

    err = arma::norm(grad, 2);  // check updated values
    if (err <= err_tol) {
        init_out_vals = x_p;
        return true;
    }

    // if (t_line < 1E-14) {
    //     printf("bfgs error: line search failed using initial values. Trying random initial values.\n");

    //     x_p.randu();
    //     t_line = line_search_mt(1.0, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, box_objfn, opt_data);
    // }

    //
    // if ||gradient(x_p)|| > tolerance, then continue

    arma::vec s = x_p - x;
    arma::vec y = grad_p - grad;

    // if (arma::norm(s,2) < 1E-14) {
    //     init_out_vals = x_p;
    //     return true;
    // }

    //
    // update W

    double W_denom_term = arma::dot(y,s);
    arma::mat W_term_1;

    if (W_denom_term > 1E-10) { // checking the curvature condition y's > 0
        W_term_1 = I_mat - s*y.t() / W_denom_term;
    
        W = W_term_1*W*W_term_1.t() + s*s.t() / W_denom_term; // update inverse Hessian approximation
    } else {
        W = 0.1*W;
    }

    grad = grad_p;

    //
    // begin loop

    int iter = 0;

    while (err > err_tol && iter < iter_max) {
        iter++;
        //
        d = - W*grad;
        line_search_mt(1.0, x_p, grad_p, d, &wolfe_cons_1, &wolfe_cons_2, box_objfn, opt_data);
        
        // err = arma::accu(arma::abs(grad_p));
        err = arma::norm(grad_p, 2);
        if (err <= err_tol) {
            break;
        }

        // if ||gradient(x_p)|| > tolerance, then continue
        s = x_p - x;
        y = grad_p - grad;

        err = arma::norm(s, 2);

        // update W
        W_denom_term = arma::dot(y,s);

        if (W_denom_term > 1E-10) { // checking the curvature condition y's > 0
            W_term_1 = I_mat - s*y.t() / W_denom_term;
        
            W = W_term_1*W*W_term_1.t() + s*s.t() / W_denom_term;
        }
        //
        x = x_p;
        grad = grad_p;
    }
    //
    if (vals_bound) {
	    x_p = inv_transform(x_p, bounds_type, lower_bounds, upper_bounds);
    }

    error_reporting(init_out_vals,x_p,opt_objfn,opt_data,success,err,err_tol,iter,iter_max,conv_failure_switch,settings_inp);
    //
    return success;
}

bool
optim::bfgs(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)> opt_objfn, void* opt_data)
{
    return bfgs_int(init_out_vals,opt_objfn,opt_data,nullptr);
}

bool
optim::bfgs(arma::vec& init_out_vals, std::function<double (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)> opt_objfn, void* opt_data, opt_settings& settings)
{
    return bfgs_int(init_out_vals,opt_objfn,opt_data,&settings);
}
