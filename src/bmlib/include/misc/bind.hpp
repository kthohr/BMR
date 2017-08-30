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
 * *bind methods; mimics R's rbind and cbind
 *
 * Keith O'Hara
 * 01/01/2012
 *
 * This version:
 * 08/14/2017
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
