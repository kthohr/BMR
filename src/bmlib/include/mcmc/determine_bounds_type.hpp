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
 * Determine the upper-lower bounds combo type
 *
 * Keith O'Hara
 * 05/01/2012
 *
 * This version:
 * 08/12/2017
 */

// note: std::isfinite is not true for: NaN, -Inf, or +Inf

inline
arma::uvec
determine_bounds_type(const bool vals_bound, const int n_vals, const arma::vec& lower_bounds, const arma::vec& upper_bounds)
{
    arma::uvec ret_vec(n_vals);

    ret_vec.fill(1); // base case: 1 - no bounds imposed

    if (vals_bound) {
        for (int i=0; i < n_vals; i++) {
            if ( std::isfinite(lower_bounds(i)) && std::isfinite(upper_bounds(i)) ) {
                // lower and upper bound imposed
                ret_vec(i) = 4;
            } else if ( std::isfinite(lower_bounds(i)) && !std::isfinite(upper_bounds(i)) ) {
                // lower bound only
                ret_vec(i) = 2;
            } else if ( !std::isfinite(lower_bounds(i)) && std::isfinite(upper_bounds(i)) ) {
                // upper bound only
                ret_vec(i) = 3;
            } // if both are infinite, this is covered by the base case
        }
    }
    //
    return ret_vec;
}
