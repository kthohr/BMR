/*################################################################################
  ##
  ##   BMLib by Keith O'Hara Copyright (C) 2011-2017
  ##   This file is part of the C++ BMLib library.
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
 * find rows with all zero elements
 *
 * Keith O'Hara
 * 01/01/2014
 *
 * This version:
 * 06/13/2017
 */

#ifndef _bmlib_zero_rows_HPP
#define _bmlib_zero_rows_HPP

arma::uvec zero_rows(const arma::mat& X);

#include "zero_rows.ipp"

#endif
