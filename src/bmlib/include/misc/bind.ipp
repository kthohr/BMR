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
 * bind methods; mimics R's rbind and cbind
 */

inline
arma::mat
cbind(const arma::mat& mat_1, const arma::mat& mat_2)
{
    return arma::join_rows(mat_1,mat_2);
}

inline
arma::mat
cbind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3)
{
    return cbind(cbind(mat_1,mat_2),mat_3);
}

inline
arma::mat
rbind(const arma::mat& mat_1, const arma::mat& mat_2)
{
    return arma::join_cols(mat_1,mat_2);
}

inline
arma::mat
rbind(const arma::mat& mat_1, const arma::mat& mat_2, const arma::mat& mat_3)
{
    return rbind(rbind(mat_1,mat_2),mat_3);
}
