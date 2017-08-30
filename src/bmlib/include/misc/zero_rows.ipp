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

inline
arma::uvec
zero_rows(const arma::mat& X)
{
    const int n = X.n_rows;
    const int k = X.n_cols;
    
    arma::mat temp_mat = arma::zeros(n,k);
    //
    temp_mat.elem(arma::find(X==0)).ones(); // find whch elements of X are = 0
    arma::colvec temp_vec = sum(temp_mat,1);
    
    arma::uvec ret_vec = arma::find(temp_vec==k); // if a row of X contains only zeros, then the corresponding row in temp_vec will = k
    //
    return ret_vec;
}
