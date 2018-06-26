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
 * A doubling algorithm to solve a discrete Lyapunov eqn
 * of the form
 *
 *            F X F' + Q = X
 *
 * The initial value for X is set to Q.
 */

//
// internal

inline
arma::mat
lyapunov_dbl_int(const arma::mat& X, const arma::mat& F, const int* max_iter_inp, const double* err_tol_inp)
{
    const int max_iter = (max_iter_inp) ? *max_iter_inp : 1000;
    const double err_tol = (err_tol_inp) ? *err_tol_inp : 1E-08;

    //

    arma::mat X_O = X; // Old
    arma::mat X_N = X; // New
    
    arma::mat dbl_mat = F;
    
    //

    int iter = 0;
    double err = 2*err_tol;

    while (err > err_tol && iter < max_iter) {
        iter++;
        
        // X_N = X_O + dbl_mat*X_O*dbl_mat.t();
        X_N += dbl_mat*X_O*dbl_mat.t();
        
        if (iter % 10 == 0) {
            err = arma::abs(X_N - X_O).max();
        }

        //

        dbl_mat *= dbl_mat;
        X_O = X_N;
    }

    //

    return X_N;
}

//
// wrappers

inline
arma::mat
lyapunov_dbl(const arma::mat& X, const arma::mat& F)
{
    return lyapunov_dbl_int(X,F,nullptr,nullptr);
}

inline
arma::mat
lyapunov_dbl(const arma::mat& X, const arma::mat& F, const int max_iter)
{
    return lyapunov_dbl_int(X,F,&max_iter,nullptr);
}

inline
arma::mat
lyapunov_dbl(const arma::mat& X, const arma::mat& F, const double err_tol)
{
    return lyapunov_dbl_int(X,F,nullptr,&err_tol);
}

inline
arma::mat
lyapunov_dbl(const arma::mat& X, const arma::mat& F, const int max_iter, const double err_tol)
{
    return lyapunov_dbl_int(X,F,&max_iter,&err_tol);
}
