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

#ifndef _bmpp_companion_form_matrix_HPP
#define _bmpp_companion_form_matrix_HPP

arma::mat companion_form_matrix(const arma::mat& mat_inp);
arma::mat companion_form_matrix(const arma::mat& mat_inp, const int c_int, const int K);

#include "companion_form_matrix.ipp"

#endif
