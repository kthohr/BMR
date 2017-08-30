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
 * Sample from a logistic distribution
 *
 * Keith O'Hara
 * 06/15/2017
 *
 * This version:
 * 07/15/2017
 */

template<typename T>
inline
T
rlogis(const T mu_par, const T sigma_par)
{
    return qlogis(runif(),mu_par,sigma_par);
}

inline
arma::mat
rlogis(const int n, const double mu_par, const double sigma_par)
{
	return rlogis(n,1,mu_par,sigma_par);
}

inline
arma::mat
rlogis(const int n, const int k, const double mu_par, const double sigma_par)
{
	arma::mat U = runif(n,k,0.0,1.0);

	return qlogis(U,mu_par,sigma_par);
}
