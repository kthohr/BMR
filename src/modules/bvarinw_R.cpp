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

RCPP_MODULE(bvarinw_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (bvarinw_R::*build_1)(const arma::mat&, bool, int) = &bvarinw_R::build_R;
    void (bvarinw_R::*build_2)(const arma::mat&, const arma::mat&, bool, int) = &bvarinw_R::build_R;

    void (bvarinw_R::*prior_3)(const arma::vec&, double, double, int) = &bvarinw_R::prior_R;

    SEXP (bvarinw_R::*forecast_1)(int, bool) = &bvarinw_R::forecast_R;
    SEXP (bvarinw_R::*forecast_2)(const arma::mat&, int, bool) = &bvarinw_R::forecast_R;
  
    // now we can declare the class
    class_<bm::bvarinw>( "bvarinw_cpp" )
        .default_constructor()

        // basic objects
        .field( "cons_term", &bm::bvarinw::cons_term )
        .field( "p", &bm::bvarinw::p )

        .field( "only_stationary_draws", &bm::bvarinw::only_stationary_draws )
        .field( "irfs_lr_restrict", &bm::bvarinw::irfs_lr_restrict )

        // read only objects
        .field_readonly( "c_int", &bm::bvarinw::c_int )
        .field_readonly( "n", &bm::bvarinw::n )
        .field_readonly( "M", &bm::bvarinw::M )
        .field_readonly( "K", &bm::bvarinw::K )
        .field_readonly( "n_ext_vars", &bm::bvarinw::n_ext_vars )

        .field_readonly( "Y", &bm::bvarinw::Y )
        .field_readonly( "X", &bm::bvarinw::X )

        .field_readonly( "alpha_hat", &bm::bvarinw::alpha_hat )
        .field_readonly( "Sigma_hat", &bm::bvarinw::Sigma_hat )

        .field_readonly( "alpha_pr_mean", &bm::bvarinw::alpha_pr_mean )
        .field_readonly( "alpha_pr_var", &bm::bvarinw::alpha_pr_var )

        .field_readonly( "Sigma_pr_scale", &bm::bvarinw::Sigma_pr_scale )
        .field_readonly( "Sigma_pr_dof", &bm::bvarinw::Sigma_pr_dof )

        .field_readonly( "alpha_pt_mean", &bm::bvarinw::alpha_pt_mean )
        .field_readonly( "alpha_pt_var", &bm::bvarinw::alpha_pt_var )

        .field_readonly( "Sigma_pt_mean", &bm::bvarinw::Sigma_pt_mean )
        .field_readonly( "Sigma_pt_dof", &bm::bvarinw::Sigma_pt_dof )

        .field_readonly( "beta_draws", &bm::bvarinw::beta_draws )
        .field_readonly( "Sigma_draws", &bm::bvarinw::Sigma_draws )
    ;

    class_<bvarinw_R>( "bvarinw" )
        .derives<bm::bvarinw>( "bvarinw_cpp" )
        .default_constructor()

        .method( "build", build_1 )
        .method( "build", build_2 )
        .method( "reset_draws", &bvarinw_R::reset_draws_R )
        .method( "prior", prior_3 )
        .method( "gibbs", &bvarinw_R::gibbs_R )
        .method( "IRF", &bvarinw_R::IRF_R )
        .method( "FEVD", &bvarinw_R::FEVD_R )
        .method( "forecast", forecast_1 )
        .method( "forecast", forecast_2 )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

void bvarinw_R::build_R(const arma::mat& data_raw, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarinw_R::build_R(const arma::mat& data_raw, const arma::mat& data_ext, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,data_ext,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarinw_R::reset_draws_R()
{
    try {
        this->reset_draws();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarinw_R::prior_R(const arma::vec& coef_prior, double Xi_beta, double Xi_Sigma, int gamma)
{
    try {
        this->prior(coef_prior,Xi_beta,Xi_Sigma,gamma);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarinw_R::gibbs_R(int n_draws, int n_burnin)
{
    try {
        this->gibbs(n_draws,n_burnin);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

SEXP bvarinw_R::IRF_R(int n_irf_periods)
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

SEXP bvarinw_R::FEVD_R(int n_periods)
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

SEXP bvarinw_R::forecast_R(int n_horizon, bool incl_shocks)
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

SEXP bvarinw_R::forecast_R(const arma::mat& Y_T, int n_horizon, bool incl_shocks)
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
