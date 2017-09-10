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
  ##   You should have received a copy of the GNU General Public License
  ##   along with BMLib. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * bind methods; mimics R's rbind and cbind
 */

#ifndef _bmlib_bind_HPP
#define _bmlib_bind_HPP

arma::mat cbind(const arma::mat& mat_1, const arma::mat& mat_2);
arma::mat cbind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3);

arma::mat rbind(const arma::mat& mat_1, const arma::mat& mat_2);
arma::mat rbind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3);

#include "bind.ipp"

#endif
