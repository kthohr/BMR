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
 *  Main gensys function, based on a discrete-time dynamic system of the form
 *
 *         Gamma0*y(t) = GammaC + Gamma1*y(t-1) + Psi*z(t) + Pi*eta(t)
 *
 *  where:         z(t) is a vector of exogenous shocks,
 *               eta(t) is a vector of 'one-step-ahead' expectational errors
 *
 *  Output is of the form:
 *
 *         y(t) = Cons + G1*y(t-1) + impact*z(t) + ywt*inv(I - fmat*inv(L))*fwt*z(t+1).
 *
 *  In a lot of cases z is iid, and so the final term drops out.
 *
 * Original code by Chris Sims
 *
 * Ported and adapted to C++ by
 * Keith O'Hara
 * 06/29/15
 *
 * This version:
 * 08/16/17
 */

#ifndef _bmpp_gensys_solver_HPP
#define _bmpp_gensys_solver_HPP

int gensys_solver(const arma::mat& Gamma0, const arma::mat& Gamma1, const arma::mat& GammaC, const arma::mat& Psi, const arma::mat& Pi,
                  arma::mat& G1, arma::mat& Cons, arma::mat& impact,
                  arma::cx_mat* fmat_out, arma::cx_mat* fwt_out, arma::cx_mat* ywt_out, arma::cx_mat* gev_out, arma::vec* eu_out, arma::mat* loose_out);

// shortened solver; only return main solver elements.
int gensys_solver(const arma::mat& Gamma0, const arma::mat& Gamma1, const arma::mat& GammaC, const arma::mat& Psi, const arma::mat& Pi,
                  arma::mat& G1, arma::mat& Cons, arma::mat& impact);


// internal
void qz_switch(const int i, arma::cx_mat& A, arma::cx_mat& B, arma::cx_mat& Q, arma::cx_mat& Z);
void qz_div(const double stake, arma::cx_mat& A, arma::cx_mat& B, arma::cx_mat& Q, arma::cx_mat& Z);

#endif
