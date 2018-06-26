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
 * Generate a sequence (of length n) from a to b.
 * Return size is n x 1 (column vector).
 */

inline
arma::vec
seqvec(const double a, const double b, const int n)
{
    if (n <= 1) {
        printf("seqvec: n needs to be > 1\n.");
        arma::vec X;

        return X;
    }

    //

    arma::vec X(n);

    if (n == 2) {
        arma::vec X(2);

        if (b >= a) {
            X(0) = a; X(1) = b;
        } else {
            X(0) = b; X(1) = a;
        }

        return X;
    }

    if (b != a) {
        double seqby = (b - a)/(n-1);
        //
        for(int ll = 0; ll < n; ll++){
            X(ll) = a + ll*seqby;
        }
    } else {
        X.fill(a);
    }

    //

    return X;
}
