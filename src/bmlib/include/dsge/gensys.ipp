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
 * Gensys
 *
 * Keith O'Hara
 * 01/01/2012
 *
 * This version:
 * 08/14/2017
 */

inline
void
gensys::build(const arma::mat& Gamma0_inp, const arma::mat& Gamma1_inp, const arma::mat& GammaC_inp, const arma::mat& Psi_inp, const arma::mat& Pi_inp)
{   
    Gamma0 = Gamma0_inp;
    Gamma1 = Gamma1_inp;
    GammaC = GammaC_inp;
    Psi    = Psi_inp;
    Pi     = Pi_inp;
}

inline
void
gensys::solve()
{
    gensys_solver(Gamma0, Gamma1, GammaC, Psi, Pi, G1, Cons, Impact);
    //
    MODEL_IS_SOLVED = true;
}

inline
void
gensys::state_space()
{
    this->state_space(F_state,G_state);

    if (F_state.n_cols > 0) {
        MODEL_HAS_STATE_FORM = true;
    }
}

inline
void
gensys::state_space(arma::mat &F_state_out, arma::mat &G_state_out)
{
    if (MODEL_IS_SOLVED) {
        arma::uvec non_expect_ind = zero_rows(Pi);
        //
        F_state_out = G1.rows(non_expect_ind);
        F_state_out = F_state_out.cols(non_expect_ind);
        
        G_state_out = Impact.rows(non_expect_ind);
    }
}

inline
arma::mat
gensys::simulate(const int sim_periods, const int burnin)
{   
    this->state_space();
    arma::mat ret_mat;

    if (shocks_cov.n_rows > 0) {
        ret_mat = dsge_simulate(F_state, G_state, shocks_cov, sim_periods, burnin);
    }

    return ret_mat;
}
