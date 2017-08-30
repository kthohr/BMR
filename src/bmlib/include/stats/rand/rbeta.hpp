/*################################################################################
  ##
  ##   Copyright (C) 2011-2017 Keith O'Hara
  ##
  ##   This file is part of the StatsLib C++ library.
  ##
  ##   StatsLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   StatsLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/* 
 * Sample from a beta distribution
 *
 * Keith O'Hara
 * 06/15/2017
 *
 * This version:
 * 07/15/2017
 */

#ifndef _statslib_rbeta_HPP
#define _statslib_rbeta_HPP

template<typename T>
T rbeta(const T a_par, const T b_par);

arma::mat rbeta(const int n, const double a_par, const double b_par);
arma::mat rbeta(const int n, const int k, const double a_par, const double b_par);

#include "rbeta.ipp"

#endif
