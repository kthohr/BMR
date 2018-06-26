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
 * Kalman filter and smoother
 */

#ifndef _bmpp_kalman_HPP
#define _bmpp_kalman_HPP

double kalman_filter(const arma::mat& data_inp, const arma::mat& F, const arma::mat& Q, const arma::mat& C, const arma::mat& H, const arma::mat& R);
double kalman_filter(const arma::mat& data_inp, const arma::mat& F, const arma::mat& Q, const arma::mat& C, const arma::mat& H, const arma::mat& R, arma::mat* state_filt_out);

#endif
