/*################################################################################
  ##
  ##   Copyright (C) 2011-2017 Keith O'Hara
  ##
  ##   This file is part of the R package BMR.
  ##
  ##   The R package BMR is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package BMR is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with BMR. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

class bvarm_R : public bm::bvarm
{
    public:
        void build_R(arma::mat data_raw, bool cons_term_inp, int p_inp);
        void build_R(arma::mat data_raw, arma::mat data_ext, bool cons_term_inp, int p_inp);

        void reset_draws_R();

        void prior_R(arma::vec coef_prior, int var_type_inp, int decay_type_inp, double HP_1_inp, double HP_2_inp, double HP_3_inp, double HP_4_inp);

        void gibbs_R(int n_draws);

        void IRF_R(int n_irf_periods);

        SEXP forecast_R(int n_horizon, bool incl_shocks);
        SEXP forecast_R(arma::mat Y_T, int n_horizon, bool incl_shocks);
};

class bvars_R : public bm::bvars
{
    public:
        void build_R(arma::mat data_raw, bool cons_term_inp, int p_inp);
        void build_R(arma::mat data_raw, arma::mat data_ext, bool cons_term_inp, int p_inp);

        void reset_draws_R();

        void prior_R(arma::vec coef_prior, double HP_1, double HP_2, arma::mat Psi_prior, double Xi_Psi, int gamma);

        void gibbs_R(int n_draws, int n_burnin);

        void IRF_R(int n_irf_periods);

        SEXP forecast_R(int n_horizon, bool incl_shocks);
        SEXP forecast_R(arma::mat Y_T, int n_horizon, bool incl_shocks);
};

class bvarw_R : public bm::bvarw
{
    public:
        void build_R(arma::mat data_raw, bool cons_term_inp, int p_inp);
        void build_R(arma::mat data_raw, arma::mat data_ext, bool cons_term_inp, int p_inp);

        void reset_draws_R();

        void prior_R(arma::vec coef_prior, double Xi_beta, double Xi_Sigma, int gamma);

        void gibbs_R(int n_draws, int n_burnin);

        void IRF_R(int n_irf_periods);

        SEXP forecast_R(int n_horizon, bool incl_shocks);
        SEXP forecast_R(arma::mat Y_T, int n_horizon, bool incl_shocks);
};

class bvartvp_R : public bm::bvartvp
{
    public:
        void build_R(arma::mat data_raw, bool cons_term_inp, int p_inp);

        void reset_draws_R();

        void prior_R(int tau_inp, double Xi_beta, double Xi_Q_inp, int gamma_Q, double Xi_Sigma, int gamma_S);

        void gibbs_R(int n_draws, int n_burnin);

        void IRF_R(int n_irf_periods, int time_period);

        // SEXP forecast_R(int n_horizon, bool incl_shocks);
        // SEXP forecast_R(arma::mat Y_T, int n_horizon, bool incl_shocks);
};

class cvar_R : public bm::cvar
{
    public:
        void build_R(arma::mat data_raw, bool cons_term_inp, int p_inp);
        void build_R(arma::mat data_raw, arma::mat data_ext, bool cons_term_inp, int p_inp);

        void reset_draws_R();

        void estim_R();

        void boot_R(int n_draws);

        void IRF_R(int n_irf_periods);

        SEXP forecast_R(int n_horizon, bool incl_shocks);
        SEXP forecast_R(arma::mat Y_T, int n_horizon, bool incl_shocks);
};
