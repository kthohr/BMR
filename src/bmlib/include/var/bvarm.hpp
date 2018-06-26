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
 * bvarm class
 */

#ifndef _bmpp_bvarm_HPP
#define _bmpp_bvarm_HPP

class bvarm
{
    public:
        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags
        int M;          // number of endogenous variables
        int K;          // number of coefficients in each
        int n_ext_vars; // number of 'external' variables

        int var_type;   // variance type (1 or 2)
        int decay_type; // decay type (1 or 2, for geometric and harmonic, resp.)

        bool only_stationary_draws = false; // discard non-stationary draws
        bool irfs_lr_restrict = false;      // impose long-run restictions on IRFs

        // 'hyper' parameters
        double HP_1;
        double HP_2;
        double HP_3;
        double HP_4;

        arma::mat Y; // Y = X beta + e
        arma::mat X;
        // arma::mat Z; // vec(Y) = Z alpha + vec(e)

        // ML-type estimates

        arma::mat alpha_hat;     // OLS estimate of alpha
        arma::mat Sigma_hat;     // covariance matrix of 'e' based on AR regressions

        // prior data

        arma::mat alpha_pr_mean; // prior mean
        arma::mat alpha_pr_var;  // prior variance

        // posterior data

        arma::vec alpha_pt_mean; // posterior mean
        arma::mat alpha_pt_var;  // posterior variance

        arma::cube beta_draws;   // posterior draws of beta

        //
        // member functions

        ~bvarm() = default;
         bvarm() = default;

        bvarm(const bvarm&) = default;
        bvarm& operator=(const bvarm&) = default;

        bvarm(bvarm&&) = default;
        bvarm& operator=(bvarm&&) = default;

        void build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp);
        void build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp);

        void reset_draws();

        void prior(const arma::vec& coef_prior, const double HP_1_inp, const double HP_2_inp, const double HP_3_inp);
        void prior(const arma::vec& coef_prior, const double HP_1_inp, const double HP_2_inp, const double HP_3_inp, const double HP_4_inp);
        void prior(const arma::vec& coef_prior, const int var_type_inp, const int decay_type_inp, const double HP_1_inp, const double HP_2_inp, const double HP_3_inp, const double HP_4_inp);

        void gibbs(const uint_t n_draws);

        arma::cube IRF(const uint_t n_irf_periods);
        arma::cube IRF(const uint_t n_irf_periods, const arma::mat& impact_mat);

        arma::cube FEVD(const uint_t n_periods);

        arma::cube forecast(const uint_t n_horizon, const bool incl_shocks);
        arma::cube forecast(const arma::mat& X_T, const uint_t n_horizon, const bool incl_shocks);

    private:
        void build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp);

        void prior_int(const arma::vec& coef_prior, const int* var_type_inp, const int* decay_type_inp, const double HP_1_inp, const double HP_2_inp, const double HP_3_inp, const double* HP_4_inp);

        arma::cube IRF_int(const uint_t n_irf_periods, const arma::mat* impact_mat);

        arma::cube forecast_int(const arma::mat* X_T_inp, const uint_t horizon, const bool incl_shocks);

        static double decay_geo(const double x, const double HP_4);
        static double decay_harm(const double x, const double HP_4);
};

//
// aux functions

// geometric decay
inline
double
bvarm::decay_geo(const double x, const double HP_4)
{
    return std::pow(HP_4,-x + 1);
}

// harmonic decay
inline
double
bvarm::decay_harm(const double x, const double HP_4)
{
    return std::pow(x, HP_4);
}

#endif
