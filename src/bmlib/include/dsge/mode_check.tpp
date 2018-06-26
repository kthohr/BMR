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
 * Modecheck
 */

template<typename T>
double
hess_objfn(const arma::vec& pars_inp, arma::vec* grad_vec, void* hess_data)
{
    mode_check_data<T>* dta = reinterpret_cast<mode_check_data<T>*>(hess_data);
    T model_obj = dta->model_obj; // need to create a copy to be thread safe

    return model_obj.log_posterior_kernel(pars_inp);
}

template<typename T>
arma::cube
mode_check(const T& model_inp, const arma::vec& mode_vals, const int grid_size)
{
    const int n_param = mode_vals.n_elem;

    mode_check_data<T> hess_data;
    hess_data.model_obj = model_inp;

    arma::mat hess_mat = optim::numerical_hessian(mode_vals,nullptr,hess_objfn<T>,&hess_data);

    arma::vec stddev_vec = arma::sqrt( - arma::diagvec( arma::inv(hess_mat) ) );

    //

    arma::mat grid_vals(grid_size,n_param), post_kernel_vals(grid_size,n_param);

    for (int j=0; j < n_param; j++) {
        grid_vals.col(j) = arma::linspace(mode_vals(j) - 2*stddev_vec(j), mode_vals(j) + 2*stddev_vec(j), grid_size);

        for (int i=0; i < grid_size; i++) {
            arma::vec pars_vec_tmp = mode_vals;
            pars_vec_tmp(j) = grid_vals(i,j);

            post_kernel_vals(i,j) = model_inp.log_posterior_kernel(pars_vec_tmp);
        }
    }

    //

    arma::cube ret_cube(grid_size,n_param,2);

    ret_cube.slice(0) = std::move(grid_vals);
    ret_cube.slice(1) = std::move(post_kernel_vals);

    return ret_cube;
}
