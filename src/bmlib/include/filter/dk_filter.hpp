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
 * Durbin-Koopman simulation smoother
 */

#ifndef _bmpp_dk_filter_HPP
#define _bmpp_dk_filter_HPP

arma::mat dk_filter(const arma::mat& Y, const arma::mat& Z, const arma::mat& Q_draw, const arma::mat& Q_chol, const arma::mat& Sigma_draw, const arma::mat& Sigma_chol, const int n_adj, const int M, const int K);

#include "dk_filter.ipp"

#endif
