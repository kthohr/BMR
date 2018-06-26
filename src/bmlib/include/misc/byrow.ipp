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
 * by-row reconstruction of a matrix ( mimics R's 'matrix(.,.,.,byrow=TRUE)' )
 */

inline
arma::mat
byrow(const arma::mat& mat_inp, const int n_rows, const int n_cols)
{
    arma::vec vec_tmp = arma::vectorise(mat_inp);
    const int n_vals = vec_tmp.n_elem;

    arma::mat ret(n_rows,n_cols);

    if (n_vals != n_rows*n_cols) { // use repmat in this case
        if (n_vals == n_rows) {
            ret = arma::repmat(vec_tmp,1,n_cols);
        } else if (n_vals == n_cols) {
            ret = arma::repmat(vec_tmp.t(),n_rows,1);
        } else {
            printf("error in byrow\n");
            return ret;
        }
    } else { // otherwise rebuild the input matrix using byrow

        int kr = 0, kc = 0;

        for (int i=0; i < n_rows*n_cols; i++) {
            ret(kr,kc) = vec_tmp(i);

            kc++;

            if (kc >= n_cols) {
                kc = 0;
                kr++;
            }
        }
    }

    return ret;
}
