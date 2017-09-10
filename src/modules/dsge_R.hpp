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

class gensys_R : public bm::gensys
{
    public:
        void build_R(const arma::mat& Gamma_0_inp, const arma::mat& Gamma_1_inp, const arma::mat& Gamma_C_inp, const arma::mat& Psi_inp, const arma::mat& Pi_inp);
        void solve_R();
        void state_space_R();
        SEXP simulate_R(int n_sim_periods, int n_burnin);
        SEXP IRF_R(int n_irf_periods);
};

class uhlig_R : public bm::uhlig
{
    public:
        void build_R(const arma::mat& A_inp, const arma::mat& B_inp, const arma::mat& C_inp, const arma::mat& D_inp, 
                     const arma::mat& F_inp, const arma::mat& G_inp, const arma::mat& H_inp, const arma::mat& J_inp, 
                     const arma::mat& K_inp, const arma::mat& L_inp, const arma::mat& M_inp, const arma::mat& N_inp);
        
        void build_pre_R(const arma::mat& A_inp, const arma::mat& B_inp, const arma::mat& C_inp, const arma::mat& D_inp);
        void build_exp_R(const arma::mat& F_inp, const arma::mat& G_inp, const arma::mat& H_inp, const arma::mat& J_inp, 
                         const arma::mat& K_inp, const arma::mat& L_inp, const arma::mat& M_inp);
        void build_exog_R(const arma::mat& N_inp);

        void solve_R();
        void state_space_R();
        SEXP simulate_R(int n_sim_periods, int n_burnin);
        SEXP IRF_R(int n_irf_periods);
};

//

class dsge_gensys_R : public bm::dsge<bm::gensys>
{
    public:
        SEXP model_fn_SEXP;
        // Rcpp::Function model_fn_Rcpp = nullptr; // crashes
        void* model_fn_Rcpp;

        arma::vec opt_initial_lb;
        arma::vec opt_initial_ub;
        arma::vec mcmc_initial_lb;
        arma::vec mcmc_initial_ub;

        void model_fn_R(const arma::vec& pars_inp, bm::gensys& lrem_obj_inp, arma::mat& shocks_cov_out, arma::mat& C_out, arma::mat& H_out, arma::mat& R_out);

        void set_model_fn(SEXP model_fn_inp);
        void eval_model(Rcpp::NumericVector pars_inp);

        void set_bounds_R(arma::vec lower_bounds_inp, arma::vec upper_bounds_inp);
        void set_prior_R(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp);

        gensys_R get_lrem_R();
        void set_lrem_R(gensys_R lrem_obj_inp);

        SEXP estim_mode_R(const arma::vec& initial_vals);
        void estim_mcmc_R(const arma::vec& initial_vals);

        void IRF_R(int n_irf_periods);
};

//

class dsgevar_gensys_R : public bm::dsgevar<bm::gensys>
{
    public:
        SEXP model_fn_SEXP;
        // Rcpp::Function model_fn_Rcpp = nullptr;
        void* model_fn_Rcpp;

        arma::vec opt_initial_lb;
        arma::vec opt_initial_ub;
        arma::vec mcmc_initial_lb;
        arma::vec mcmc_initial_ub;

        // member functions

        void build_R(const arma::mat& data_raw, bool cons_term_inp, int p_inp, double lambda_inp);

        arma::mat get_dsge_draws();

        void model_fn_R(const arma::vec& pars_inp, bm::gensys& lrem_obj_inp, arma::mat& shocks_cov_out, arma::mat& C_out, arma::mat& H_out, arma::mat& R_out);

        void set_model_fn(SEXP model_fn_inp);
        void eval_model(Rcpp::NumericVector pars_inp);

        void set_bounds_R(arma::vec lower_bounds_inp, arma::vec upper_bounds_inp);
        void set_prior_R(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp);

        gensys_R get_lrem_R();
        void set_lrem_R(gensys_R lrem_obj_inp);

        dsge_gensys_R get_dsge_R();
        void set_dsge_R(dsge_gensys_R dsge_obj_inp);

        SEXP estim_mode_R(const arma::vec& initial_vals);
        void estim_mcmc_R(const arma::vec& initial_vals);

        void IRF_R(int n_irf_periods);
};
