/*################################################################################
  ##
  ##   MacroLib by Keith O'Hara Copyright (C) 2011-2016
  ##   This file is part of the C++ MacroLib library.
  ##
  ##   MacroLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   MacroLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/*
 * Sample from a mutivariate normal distribution
 * mu is K x 1 vector of zeros; Sigma is K x K
 * returns a K x 1 vector
 */

arma::vec rmvnorm(arma::mat Sigma, bool pre_chol)
{
	int K = Sigma.n_cols;
	arma::mat A;

	if (pre_chol) {
		A = Sigma;
	} else {
		A = arma::trans(arma::chol(Sigma));
	}

	arma::vec ret = A * arma::randn(K);
	//
	return ret;
}

/*
 * Sample from a mutivariate normal distribution
 * mu is K x 1; Sigma is K x K
 * returns a K x 1 vector
 */

arma::vec rmvnorm(arma::mat mu, arma::mat Sigma, bool pre_chol)
{
	int K = Sigma.n_cols;
	arma::mat A;

	if (pre_chol) {
		A = Sigma;
	} else {
		A = arma::trans(arma::chol(Sigma));
	}

	arma::vec ret = mu + A * arma::randn(K);
	//
	return ret;
}

/*
 * Sample from a mutivariate normal distribution
 * mu is K x 1 vector of zeros; Sigma is K x K
 * returns an n x K vector
 */

arma::mat rmvnorm(int n, arma::mat Sigma, bool pre_chol)
{
	int K = Sigma.n_cols;
	arma::mat A;

	if (pre_chol) {
		A = Sigma;
	} else {
		A = arma::trans(arma::chol(Sigma));
	}
	
	arma::mat ret = arma::trans(A * arma::randn(K,n));
	//
	return ret;
}

/*
 * Sample from a mutivariate normal distribution
 * mu is K x 1; Sigma is K x K
 * returns an n x K vector
 */

arma::mat rmvnorm(int n, arma::mat mu, arma::mat Sigma, bool pre_chol)
{
	int K = Sigma.n_cols;
	arma::mat A;

	if (pre_chol) {
		A = Sigma;
	} else {
		A = arma::trans(arma::chol(Sigma));
	}
	
	arma::mat ret = arma::repmat(mu.t(),n,1) + arma::trans(A * arma::randn(K,n));
	//
	return ret;
}
