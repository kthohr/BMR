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
 * cube to matrix
 */

inline
arma::mat
cube_to_mat(const arma::cube& cube_inp)
{
    return cube_to_mat(cube_inp,false);
}

inline
arma::mat
cube_to_mat(const arma::cube& cube_inp, const bool vec_op)
{
    const int cube_rows   = cube_inp.n_rows;
    const int cube_cols   = cube_inp.n_cols;
    const int cube_slices = cube_inp.n_slices;

    //

    int mat_rows = 0;
    int mat_cols = 0;

    if (vec_op) {
        mat_rows = cube_slices;
        mat_cols = cube_rows*cube_cols;
    } else {
        mat_rows = cube_rows*cube_slices;
        mat_cols = cube_cols;
    }

    //

    arma::mat ret_mat(mat_rows,mat_cols);

    for (int i=0; i < cube_slices; i++) {
        if (vec_op) {
            ret_mat.row(i) = arma::trans(arma::vectorise(cube_inp.slice(i)));
        } else {
            ret_mat.rows(i*cube_rows,(i+1)*cube_rows - 1) = cube_inp.slice(i);
        }
    }

    //

    return ret_mat;
}
