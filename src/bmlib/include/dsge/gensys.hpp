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

#ifndef _bmlib_gensys_HPP
#define _bmlib_gensys_HPP

class gensys
{
    public:
        // solution input matrices
        arma::mat Gamma0;
        arma::mat Gamma1;
        arma::mat GammaC;
        arma::mat Psi;
        arma::mat Pi;
         
        // solution matrices
        arma::mat G1;
        arma::mat Impact;
        arma::mat Cons;
         
        // state-space form
        arma::mat F_state;
        arma::mat G_state;
         
        // covariance matrix of shocks
        arma::mat shocks_cov;
         
        // member functions
        void build(const arma::mat& Gamma0_inp, const arma::mat& Gamma1_inp, const arma::mat& GammaC_inp, const arma::mat& Psi_inp, const arma::mat& Pi_inp);

        void solve();

        arma::mat simulate(const int sim_periods, const int burnin);

        void state_space();
        void state_space(arma::mat& F_state_out, arma::mat& G_state_out);
         
    protected:
        bool MODEL_IS_SOLVED      = false;
        bool MODEL_HAS_STATE_FORM = false;
};

#include "gensys.ipp"

#endif
