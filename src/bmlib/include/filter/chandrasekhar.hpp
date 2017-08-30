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
 * Chandrasekhar recursions
 *
 * Keith O'Hara
 * 01/01/2012
 *
 * This version:
 * 08/14/2017
 */

#ifndef _bmlib_chandrasekhar_HPP
#define _bmlib_chandrasekhar_HPP

double chand_recur(const arma::mat& data_inp, const arma::mat& F, const arma::mat& Q, const arma::mat& C, const arma::mat& H, const arma::mat& R);

#endif
