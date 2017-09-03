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
 * uhlig class
 */

inline
void uhlig::solve()
{
    uhlig_solver(A,B,C,D,F,G,H,J,K,L,M,N,P_sol,Q_sol,R_sol,S_sol);
}

inline
void
uhlig::state_space()
{
    this->state_space(F_state,G_state);
}

inline
void 
uhlig::state_space(arma::mat &F_state_out, arma::mat &G_state_out)
{
    const int dim_k = N.n_rows;
    const int dim_m = P_sol.n_rows;
    const int dim_n = R_sol.n_rows;

    F_state_out = rbind( cbind(rbind(P_sol,R_sol),arma::zeros(dim_m+dim_n,dim_n),rbind(Q_sol,S_sol)*N), cbind(arma::zeros(dim_k,dim_m + dim_n), N) );
    G_state_out = rbind(Q_sol,S_sol,arma::eye(dim_k,dim_k));
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

inline
arma::cube
uhlig::IRF(const int n_irf_periods)
{
    this->state_space();
    arma::cube ret_cube;

    if (shocks_cov.n_rows > 0) {
        ret_cube = dsge_irf(F_state, G_state, shocks_cov, n_irf_periods);
    }

    return ret_cube;
}
