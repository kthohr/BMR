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
 * uhlig class
 */

#ifndef _bmpp_uhlig_HPP
#define _bmpp_uhlig_HPP

class uhlig
{
    public:
        // input matrices
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
        arma::mat P_sol;
        arma::mat Q_sol;
        arma::mat R_sol;
        arma::mat S_sol;

        // state-space form
        arma::mat F_state;
        arma::mat G_state;

        // covariance matrix of shocks
        arma::mat shocks_cov;

        //
        // member functions

        ~uhlig() = default;
         uhlig() = default;

        uhlig(const uhlig&) = default;
        uhlig& operator=(const uhlig&) = default;

        uhlig(uhlig&&) = default;
        uhlig& operator=(uhlig&&) = default;

        void build(const arma::mat& A_inp, const arma::mat& B_inp, const arma::mat& C_inp, const arma::mat& D_inp, 
                   const arma::mat& F_inp, const arma::mat& G_inp, const arma::mat& H_inp, const arma::mat& J_inp, 
                   const arma::mat& K_inp, const arma::mat& L_inp, const arma::mat& M_inp, const arma::mat& N_inp);

        void build_pre(const arma::mat& A_inp, const arma::mat& B_inp, const arma::mat& C_inp, const arma::mat& D_inp);
        void build_exp(const arma::mat& F_inp, const arma::mat& G_inp, const arma::mat& H_inp, const arma::mat& J_inp, 
                       const arma::mat& K_inp, const arma::mat& L_inp, const arma::mat& M_inp);
        void build_exog(const arma::mat& N_inp);

        void solve();

        void state_space();
        void state_space(arma::mat& F_state, arma::mat& G_state);

        arma::mat simulate(const int sim_periods, const int burnin);

        arma::cube IRF(const int n_irf_periods);
};

#include "uhlig.ipp"

#endif
