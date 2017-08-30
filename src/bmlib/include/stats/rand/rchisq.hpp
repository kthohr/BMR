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
 * Sample from a Chi-Squared distribution
 *
 * Keith O'Hara
 * 06/01/2015
 *
 * This version:
 * 07/15/2017
 */

#ifndef _statslib_rchisq_HPP
#define _statslib_rchisq_HPP

double rchisq(const int dof_par);

arma::mat rchisq(const int n, const int dof_par);
arma::mat rchisq(const int n, const int k, const int dof_par);

#include "rchisq.ipp"

#endif
