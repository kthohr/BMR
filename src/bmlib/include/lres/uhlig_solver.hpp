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
 * Uhlig's method for solving linear rational  expectations
 * models using a generalized eigen decomposition.
 *
 * Written by
 * Keith O'Hara
 * 07/01/12
 *
 * This version:
 * 01/03/17
 *
 * Requires the 'Armadillo' linear algebra library.
 */

int uhlig_solver(const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& D,
                 const arma::mat& F, const arma::mat& G, const arma::mat& H, const arma::mat& J, const arma::mat& K, const arma::mat& L, const arma::mat& M, const arma::mat& N,
                 const arma::vec* which_eig, arma::mat& P, arma::mat& Q, arma::mat& R, arma::mat& S, arma::cx_vec* eigen_vals_out, arma::cx_mat* eigen_vecs_out);

// without which_eig or returning the eigenvalues and eigenvectors
int uhlig_solver(const arma::mat& A, const arma::mat& B, const arma::mat& C, const arma::mat& D,
                 const arma::mat& F, const arma::mat& G, const arma::mat& H, const arma::mat& J, const arma::mat& K, const arma::mat& L, const arma::mat& M, const arma::mat& N,
                 arma::mat& P, arma::mat& Q, arma::mat& R, arma::mat& S);
