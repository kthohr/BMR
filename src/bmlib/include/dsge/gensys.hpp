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

#ifndef _bmpp_gensys_HPP
#define _bmpp_gensys_HPP

class gensys
{
    public:
        // input matrices: Gamma0*X_t = GammaC + Gamma1*X_{t-1} + Psi*z_t + Pi*eta_t
        arma::mat Gamma_0;
        arma::mat Gamma_1;
        arma::mat Gamma_C;
        arma::mat Psi;
        arma::mat Pi;

        // solution matrices: X_t = C + G*X_{t-1} + imp*z_t
        arma::mat G_sol;
        arma::mat impact_sol;
        arma::mat cons_sol;

        // state-space form
        arma::mat F_state;
        arma::mat G_state;

        // covariance matrix of shocks
        arma::mat shocks_cov;

        //
        // member functions

        ~gensys() = default;
         gensys() = default;

        gensys(const gensys&) = default;
        gensys& operator=(const gensys&) = default;

        gensys(gensys&&) = default;
        gensys& operator=(gensys&&) = default;

        void build(const arma::mat& Gamma_0_inp, const arma::mat& Gamma_1_inp, const arma::mat& Gamma_C_inp, const arma::mat& Psi_inp, const arma::mat& Pi_inp);

        void solve();

        void state_space();
        void state_space(arma::mat& F_state_out, arma::mat& G_state_out);

        arma::mat simulate(const int n_sim_periods, const int n_burnin);

        arma::cube IRF(const int n_irf_periods);
};

#include "gensys.ipp"

#endif
