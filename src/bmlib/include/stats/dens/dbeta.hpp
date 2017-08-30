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
 * pdf of the beta distribution
 *
 * Keith O'Hara
 * 01/02/16
 *
 * This version:
 * 07/06/2017
 */

#ifndef _statslib_dbeta_HPP
#define _statslib_dbeta_HPP

// single input
template<typename T>
statslib_constexpr T dbeta(const T x, const T a_par, const T b_par, const bool log_form);

statslib_constexpr double dbeta(const double x);
statslib_constexpr double dbeta(const double x, const bool log_form);
statslib_constexpr double dbeta(const double x, const double a_par, const double b_par);

// matrix/vector input
arma::mat dbeta_int(const arma::mat& x, const double* a_par_inp, const double* b_par_inp, const bool log_form);

arma::mat dbeta(const arma::mat& x);
arma::mat dbeta(const arma::mat& x, const bool log_form);
arma::mat dbeta(const arma::mat& x, const double a_par, const double b_par);
arma::mat dbeta(const arma::mat& x, const double a_par, const double b_par, const bool log_form);

#include "dbeta.ipp"

#endif
