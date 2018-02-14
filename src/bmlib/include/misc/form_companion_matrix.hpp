/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
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
  ##   You should have received a copy of the GNU General Public License
  ##   along with BMLib. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * build a companion form matrix
 */

inline
arma::mat
form_companion_matrix(const arma::mat& A, const int K, const int p)
{
    arma::mat A_comp;

    if (p > 1) {
        A_comp = arma::join_cols(arma::trans(A), arma::join_rows(arma::eye(K*(p-1),K*(p-1)), arma::zeros(K*(p-1), K)));
    } else {
        A_comp = A;
    }

    return A_comp;
}
