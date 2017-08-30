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
 * Simulate data from a DSGE Model
 *
 * Keith O'Hara
 * 01/01/2012
 *
 * This version:
 * 08/15/2017
 */

inline
arma::mat
dsge_simulate(const arma::mat& F_state, const arma::mat& G_state, const arma::mat& shocks_cov, const int sim_periods, const int burnin)
{
    // const int n_shocks = G_state.n_cols;

    arma::mat dsge_sim_mat(sim_periods+burnin,F_state.n_cols);

    arma::mat shocks = stats::rmvnorm(sim_periods + burnin, shocks_cov);
    //
    dsge_sim_mat.row(0) = shocks.row(0) * G_state.t();

    for (int i = 1; i < (sim_periods + burnin); i++) {
        dsge_sim_mat.row(i) = dsge_sim_mat.row(i-1) * F_state.t() + shocks.row(i) * G_state.t();
    }

    dsge_sim_mat.shed_rows(0,burnin-1);
    //
    return dsge_sim_mat;
}
