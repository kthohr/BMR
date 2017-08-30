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

#ifndef _bmlib_uhlig_HPP
#define _bmlib_uhlig_HPP

class uhlig
{
    public:
        // solution input matrices
        arma::mat A;
        arma::mat B;
        arma::mat C;
        arma::mat D;
        
        arma::mat F;
        arma::mat G;
        arma::mat H;
        
        arma::mat J;
        arma::mat K;
        arma::mat L;
        arma::mat M;
        arma::mat N;
        
        // solution matrices
        arma::mat P;
        arma::mat Q;
        arma::mat R;
        arma::mat S;
        
        // state-space form
        arma::mat F_state;
        arma::mat G_state;
        
        // covariance matrix of shocks
        arma::mat shocks_cov;
        
        // member functions
        void build_block_pre(const arma::mat* A_inp, const arma::mat* B_inp, const arma::mat* C_inp, const arma::mat* D_inp);
        void build_block_exp(const arma::mat* F_inp, const arma::mat* G_inp, const arma::mat* H_inp, const arma::mat* J_inp, const arma::mat* K_inp, const arma::mat* L_inp, const arma::mat* M_inp);
        void build_block_exog(const arma::mat* N_inp);

        void solve();

        arma::mat simulate(const int sim_periods, const int burnin);

        void state_space();
        void state_space(arma::mat& F_state, arma::mat& G_state);
        
    private:
        bool MODEL_IS_SOLVED      = false;
        bool MODEL_HAS_STATE_FORM = false;
};

#include "uhlig.ipp"

#endif
