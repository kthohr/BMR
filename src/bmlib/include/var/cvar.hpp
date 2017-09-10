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
  ##   You should have received a copy of the GNU General Public License
  ##   along with BMLib. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * cvar class
 */

class cvar
{
    public:
        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags
        int M;          // number of endogenous variables
        int K;          // number of coefficients in each
        int n_ext_vars; // number of 'external' variables

        arma::mat Y;    // Y = X beta + e
        arma::mat X;

        arma::mat beta_hat;      // OLS estimate of beta
        arma::mat Sigma_hat;     // OLS-based estimation of covariance matrix of 'e'

        arma::cube beta_draws;   // bootstrap draws of beta
        arma::cube Sigma_draws;  // bootstrap draws of beta

        arma::cube irfs;         // irfs based on the bootstrap draws

        //
        // member functions

        ~cvar() = default;
         cvar() = default;

        cvar(const cvar&) = default;
        cvar& operator=(const cvar&) = default;

        cvar(cvar&&) = default;
        cvar& operator=(cvar&&) = default;

        void build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp);
        void build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp);

        void reset_draws();

        void estim();

        void boot(const int n_draws);

        void IRF(const int n_irf_periods);

        arma::cube forecast(const int horizon, const bool incl_shocks);
        arma::cube forecast(const arma::mat& Y_T, const int horizon, const bool incl_shocks);

    private:
        void build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp);
        arma::cube forecast_int(const arma::mat* Y_T_inp, const int horizon, const bool incl_shocks);
};
