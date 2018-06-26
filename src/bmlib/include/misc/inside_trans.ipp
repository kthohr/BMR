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
 * transpose square blocks
 */

inline
arma::mat
inside_trans(const arma::mat& mat_inp, const bool skip_first_row)
{
    const int c_int = !!skip_first_row;

    const int K = mat_inp.n_rows;
    const int M = mat_inp.n_cols;

    const int p = (K-c_int) / M;

    //

    arma::mat ret_mat = mat_inp;

    for (int i=0; i < p; i++) {
        ret_mat.rows(c_int+i*M,(i+1)*M) = ret_mat.rows(c_int+i*M,(i+1)*M).t();
    }

    return ret_mat;
}
