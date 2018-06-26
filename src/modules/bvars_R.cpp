/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
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

#include "bmpp.hpp"
#include "vars_R.hpp"

RCPP_MODULE(bvars_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (bvars_R::*build_1)(const arma::mat&, bool, int) = &bvars_R::build_R;
    void (bvars_R::*build_2)(const arma::mat&, const arma::mat&, bool, int) = &bvars_R::build_R;

    void (bvars_R::*prior_3)(const arma::vec&, double, double, const arma::mat&, double, int, bool) = &bvars_R::prior_R;

    SEXP (bvars_R::*forecast_1)(int, bool) = &bvars_R::forecast_R;
    SEXP (bvars_R::*forecast_2)(const arma::mat&, int, bool) = &bvars_R::forecast_R;
  
    // now we can declare the class
    class_<bm::bvars>( "bvars_cpp" )
        .default_constructor()

        // basic objects
        .field( "cons_term", &bm::bvars::cons_term )
        .field( "p", &bm::bvars::p )

        .field( "only_stationary_draws", &bm::bvars::only_stationary_draws )
        .field( "irfs_lr_restrict", &bm::bvars::irfs_lr_restrict )

        // read only objects
        .field_readonly( "c_int", &bm::bvars::c_int )
        .field_readonly( "n", &bm::bvars::n )
        .field_readonly( "M", &bm::bvars::M )
        .field_readonly( "K", &bm::bvars::K )
        .field_readonly( "n_ext_vars", &bm::bvars::n_ext_vars )

        .field_readonly( "Y", &bm::bvars::Y )
        .field_readonly( "X", &bm::bvars::X )

        .field_readonly( "psi_hat", &bm::bvars::psi_hat )
        .field_readonly( "alpha_hat", &bm::bvars::alpha_hat )
        .field_readonly( "Sigma_hat", &bm::bvars::Sigma_hat )

        .field_readonly( "psi_pr_mean", &bm::bvars::psi_pr_mean )
        .field_readonly( "psi_pr_var", &bm::bvars::psi_pr_var )

        .field_readonly( "alpha_pr_mean", &bm::bvars::alpha_pr_mean )
        .field_readonly( "alpha_pr_var", &bm::bvars::alpha_pr_var )

        .field_readonly( "Sigma_pr_scale", &bm::bvars::Sigma_pr_scale )
        .field_readonly( "Sigma_pr_dof", &bm::bvars::Sigma_pr_dof )

        .field_readonly( "psi_pt_mean", &bm::bvars::psi_pt_mean )
        .field_readonly( "psi_pt_var", &bm::bvars::psi_pt_var )

        .field_readonly( "alpha_pt_mean", &bm::bvars::alpha_pt_mean )
        .field_readonly( "alpha_pt_var", &bm::bvars::alpha_pt_var )

        .field_readonly( "Sigma_pt_mean", &bm::bvars::Sigma_pt_mean )
        .field_readonly( "Sigma_pt_dof", &bm::bvars::Sigma_pt_dof )

        .field_readonly( "Psi_draws", &bm::bvars::Psi_draws )
        .field_readonly( "beta_draws", &bm::bvars::beta_draws )
        .field_readonly( "Sigma_draws", &bm::bvars::Sigma_draws )
    ;

    class_<bvars_R>( "bvars" )
        .derives<bm::bvars>( "bvars_cpp" )
        .default_constructor()

        .method( "build", build_1 )
        .method( "build", build_2 )
        .method( "reset_draws", &bvars_R::reset_draws_R )
        .method( "prior", prior_3 )
        .method( "gibbs", &bvars_R::gibbs_R )
        .method( "IRF", &bvars_R::IRF_R )
        .method( "FEVD", &bvars_R::FEVD_R )
        .method( "forecast", forecast_1 )
        .method( "forecast", forecast_2 )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

void bvars_R::build_R(const arma::mat& data_raw, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvars_R::build_R(const arma::mat& data_raw, const arma::mat& data_ext, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,data_ext,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvars_R::reset_draws_R()
{
    try {
        this->reset_draws();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvars_R::prior_R(const arma::vec& coef_prior, double HP_1, double HP_2, 
                      const arma::mat& Psi_prior, double Xi_psi, int gamma,
                      bool full_cov_prior)
{
    try {
        this->prior(coef_prior,HP_1,HP_2,Psi_prior,Xi_psi,gamma,full_cov_prior);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvars_R::gibbs_R(int n_draws, int n_burnin)
{
    try {
        this->gibbs(n_draws,n_burnin);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

SEXP bvars_R::IRF_R(int n_irf_periods)
{
    try {
        arma::cube irf_vals = this->IRF(n_irf_periods);

        return Rcpp::List::create(Rcpp::Named("irf_vals") = irf_vals);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP bvars_R::FEVD_R(int n_periods)
{
    try {
        arma::cube fevd_vals = this->FEVD(n_periods);

        return Rcpp::List::create(Rcpp::Named("fevd_vals") = fevd_vals);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP bvars_R::forecast_R(int n_horizon, bool incl_shocks)
{
    try {
        arma::cube fcast_res = this->forecast(n_horizon,incl_shocks);

        return Rcpp::List::create(Rcpp::Named("forecast_vals") = fcast_res);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP bvars_R::forecast_R(const arma::mat& Y_T, int n_horizon, bool incl_shocks)
{
    try {
        arma::cube fcast_res = this->forecast(Y_T,n_horizon,incl_shocks);

        return Rcpp::List::create(Rcpp::Named("forecast_vals") = fcast_res);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
