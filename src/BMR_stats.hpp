/*################################################################################
  ##
  ##   R package BMR by Keith O'Hara Copyright (C) 2011-2016
  ##   This file is part of the R package BMR.
  ##
  ##   The R package BMR is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package BMR is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

#ifndef _BMR_stats_HPP
#define _BMR_stats_HPP

#include <RcppArmadillo.h>

double rchisq(int k);
arma::colvec rchisq(int n, int k);

arma::mat rinvwish(arma::mat Psi, int nu);

arma::vec rmvnorm(arma::mat Sigma, bool pre_chol);
arma::vec rmvnorm(arma::mat mu, arma::mat Sigma, bool pre_chol);
arma::mat rmvnorm(int n, arma::mat Sigma, bool pre_chol);
arma::mat rmvnorm(int n, arma::mat mu, arma::mat Sigma, bool pre_chol);

#endif