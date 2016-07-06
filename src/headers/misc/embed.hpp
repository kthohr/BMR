/*################################################################################
  ##
  ##   MacroLib by Keith O'Hara Copyright (C) 2011-2016
  ##   This file is part of the C++ MacroLib library.
  ##
  ##   MacroLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   MacroLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/*
 * embed
 */

arma::mat embed(arma::mat X, int p)
{
    int i;
    int n = X.n_rows;
    int k = X.n_cols;

    arma::mat ret(n-p,k*(p+1));
    //
    for (i=0; i<(p+1); i++) {
        ret.cols(i*k,(i+1)*k-1) = X.rows(p-i,n-i-1);
    }
    //
    return ret;
}
