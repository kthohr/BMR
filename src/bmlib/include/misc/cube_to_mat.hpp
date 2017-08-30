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
 * cube to matrix
 *
 * Keith O'Hara
 * 01/01/2012
 *
 * This version:
 * 08/28/2017
 */

#ifndef _cube_to_mat_HPP
#define _cube_to_mat_HPP

arma::mat cube_to_mat(const arma::cube& cube_inp);
arma::mat cube_to_mat(const arma::cube& cube_inp, const bool vec_op);

#include "cube_to_mat.ipp"

#endif
