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
 * gensys class
 */

inline
void
gensys::build(const arma::mat& Gamma_0_inp, const arma::mat& Gamma_1_inp, const arma::mat& Gamma_C_inp, const arma::mat& Psi_inp, const arma::mat& Pi_inp)
{
    Gamma_0 = Gamma_0_inp;
    Gamma_1 = Gamma_1_inp;
    Gamma_C = Gamma_C_inp;
    Psi = Psi_inp;
    Pi  = Pi_inp;
}

inline
void
gensys::solve()
{
    gensys_solver(Gamma_0, Gamma_1, Gamma_C, Psi, Pi, G_sol, cons_sol, impact_sol);
}

inline
void
gensys::state_space()
{
    this->state_space(F_state,G_state);
}

inline
void
gensys::state_space(arma::mat &F_state_out, arma::mat &G_state_out)
{
    if (impact_sol.n_elem > 0) {
        arma::uvec non_expect_ind = zero_rows(Pi);
        //
        F_state_out = G_sol.rows(non_expect_ind);
        F_state_out = F_state_out.cols(non_expect_ind);

        G_state_out = impact_sol.rows(non_expect_ind);
    }
}

inline
arma::mat
gensys::simulate(const int n_sim_periods, const int n_burnin)
{
    this->state_space();
    arma::mat ret_mat;

    if (shocks_cov.n_rows > 0) {
        ret_mat = dsge_simulate(F_state, G_state, shocks_cov, n_sim_periods, n_burnin);
    } else {
        printf("dsge simulate: error: shocks_cov not found\n");
    }

    return ret_mat;
}

inline
arma::cube
gensys::IRF(const int n_irf_periods)
{
    this->state_space();
    arma::cube ret_cube;

    if (shocks_cov.n_rows > 0) {
        ret_cube = dsge_irf(F_state, G_state, shocks_cov, n_irf_periods);
    } else {
        printf("dsge irf: error: shocks_cov not found\n");
    }

    return ret_cube;
}
