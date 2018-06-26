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
 * DSGE IRF
 */

inline
arma::cube
dsge_irf(const arma::mat& F_state, const arma::mat& G_state, const arma::mat& shocks_cov, const int n_irf_periods)
{
    const int n_shocks = G_state.n_cols;
    const int n_vars = F_state.n_cols;

    arma::cube dsge_irf_out(n_irf_periods,n_vars,n_shocks);

    //

    for (int j=0; j < n_shocks; j++) {
        
        arma::rowvec shock_vec = arma::zeros(1,n_shocks);
        shock_vec(j) = std::sqrt(shocks_cov(j,j)); // 1 std. dev. shock

        arma::mat irf_mat = arma::zeros(n_irf_periods,n_vars);

        irf_mat.row(0) = shock_vec*G_state.t();

        if (n_irf_periods > 1) {
            for (int i=1; i < n_irf_periods; i++) {
                irf_mat.row(i) = irf_mat.row(i-1)*F_state.t();
            }
        }

        dsge_irf_out.slice(j) = irf_mat;
    }

    //

    return dsge_irf_out;
}
