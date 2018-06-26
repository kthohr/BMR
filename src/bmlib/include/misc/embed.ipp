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
 * Embed a matrix
 */

inline 
arma::mat 
embed(const arma::mat& X, const int p)
{
    const int n = X.n_rows;
    const int k = X.n_cols;

    //

    arma::mat ret(n-p,k*(p+1));
    
    for (int i=0; i < (p+1); i++) {
        ret.cols(i*k,(i+1)*k-1) = X.rows(p-i,n-i-1);
    }

    //

    return ret;
}

inline
void
embed(arma::mat &X, arma::mat &Y, const int p)
{
	int n = X.n_rows;
	int K = X.n_cols;

    //
	
	arma::mat Z(n-p,K*(p+1));

	for (int j=0; j <= p; j++)
    {
		Z.cols(j*K,(j+1)*K-1) = X.rows(p-j,n-1-j);
	}

	//
    
	Y = Z.cols(0,K-1);
	X = Z.cols(K,K*(p+1)-1);
}
