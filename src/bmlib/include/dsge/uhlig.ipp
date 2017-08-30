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
 * Uhlig
 *
 * Keith O'Hara
 * 01/01/2012
 *
 * This version:
 * 08/14/2017
 */

inline
void uhlig::solve()
{
    uhlig_solver(A,B,C,D,F,G,H,J,K,L,M,N,P,Q,R,S);
    //
    MODEL_IS_SOLVED = true;
}

inline
void
uhlig::state_space()
{
    this->state_space(F_state,G_state);

    if (F_state.n_cols > 0) {
        MODEL_HAS_STATE_FORM = true;
    }
}

inline
void 
uhlig::state_space(arma::mat &F_state_out, arma::mat &G_state_out)
{
    if (MODEL_IS_SOLVED) {
        const int dim_k = N.n_rows;
        const int dim_m = P.n_rows;
        const int dim_n = R.n_rows;

        F_state_out = rbind( cbind(rbind(P,R),arma::zeros(dim_m+dim_n,dim_n),rbind(Q,S)*N), cbind(arma::zeros(dim_k,dim_m + dim_n), N) );
        G_state_out = rbind(Q,S,arma::eye(dim_k,dim_k));
    }
}

inline
arma::mat
uhlig::simulate(const int sim_periods, const int burnin)
{   
    this->state_space();
    arma::mat ret_mat;

    if (shocks_cov.n_rows > 0) {
        ret_mat = dsge_simulate(F_state, G_state, shocks_cov, sim_periods, burnin);
    }

    return ret_mat;
}
