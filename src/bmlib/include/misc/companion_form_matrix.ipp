/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the BM++ C++ library.
  ##
  ##   BM++ is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   BM++ is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with BM++. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * build a companion form matrix: Y = A X + e
 */

inline
arma::mat
companion_form_matrix(const arma::mat& mat_inp)
{
    const uint_t M = mat_inp.n_cols;
    const uint_t p = mat_inp.n_rows / M; // note that 3/2, etc cases, will default to 'floor'

    arma::mat A_comp;

    if (p > 1) {
        const uint_t Mpm1 = M*(p-1);
        A_comp = arma::join_cols(mat_inp.t(), arma::join_rows(arma::eye(Mpm1,Mpm1), arma::zeros(Mpm1,M)));
    } else {
        A_comp = mat_inp.t();
    }

    return A_comp;
}

inline
arma::mat
companion_form_matrix(const arma::mat& mat_inp, const int c_int, const int K)
{
    return companion_form_matrix(mat_inp.rows(c_int,K-1));
}
