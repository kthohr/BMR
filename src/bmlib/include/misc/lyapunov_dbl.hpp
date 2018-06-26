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
 * A doubling algorithm to solve a discrete Lyapunov eqn
 * of the form
 *
 *            F X F' + Q = X
 *
 * The initial value for X is set to Q.
 */

 #ifndef _bmpp_lyapunov_dbl_HPP
 #define _bmpp_lyapunov_dbl_HPP

// internal
arma::mat lyapunov_dbl_int(const arma::mat& X, const arma::mat& F, const int* max_iter_inp, const double* err_tol_inp);

// wrappers
arma::mat lyapunov_dbl(const arma::mat& X, const arma::mat& F);
arma::mat lyapunov_dbl(const arma::mat& X, const arma::mat& F, const int max_iter);
arma::mat lyapunov_dbl(const arma::mat& X, const arma::mat& F, const double err_tol);
arma::mat lyapunov_dbl(const arma::mat& X, const arma::mat& F, const int max_iter, const double err_tol);

#include "lyapunov_dbl.ipp"

#endif
