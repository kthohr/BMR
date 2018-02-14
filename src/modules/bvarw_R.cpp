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

#include "bmlib.hpp"
#include "vars_R.hpp"

RCPP_MODULE(bvarw_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (bvarw_R::*build_1)(const arma::mat&, bool, int) = &bvarw_R::build_R;
    void (bvarw_R::*build_2)(const arma::mat&, const arma::mat&, bool, int) = &bvarw_R::build_R;

    void (bvarw_R::*prior_3)(const arma::vec&, double, double, int) = &bvarw_R::prior_R;

    SEXP (bvarw_R::*forecast_1)(int, bool) = &bvarw_R::forecast_R;
    SEXP (bvarw_R::*forecast_2)(const arma::mat&, int, bool) = &bvarw_R::forecast_R;
  
    // now we can declare the class
    class_<bm::bvarw>( "bvarw_cpp" )
        .default_constructor()

        // basic objects
        .field( "cons_term", &bm::bvarw::cons_term )
        .field( "p", &bm::bvarw::p )

        .field( "only_stationary_draws", &bm::bvarw::only_stationary_draws )
        .field( "irfs_lr_restrict", &bm::bvarw::irfs_lr_restrict )

        // read only objects
        .field_readonly( "c_int", &bm::bvarw::c_int )
        .field_readonly( "n", &bm::bvarw::n )
        .field_readonly( "M", &bm::bvarw::M )
        .field_readonly( "K", &bm::bvarw::K )
        .field_readonly( "n_ext_vars", &bm::bvarw::n_ext_vars )

        .field_readonly( "Y", &bm::bvarw::Y )
        .field_readonly( "X", &bm::bvarw::X )

        .field_readonly( "alpha_hat", &bm::bvarw::alpha_hat )
        .field_readonly( "Sigma_hat", &bm::bvarw::Sigma_hat )

        .field_readonly( "alpha_pr_mean", &bm::bvarw::alpha_pr_mean )
        .field_readonly( "alpha_pr_var", &bm::bvarw::alpha_pr_var )

        .field_readonly( "Sigma_pr_scale", &bm::bvarw::Sigma_pr_scale )
        .field_readonly( "Sigma_pr_dof", &bm::bvarw::Sigma_pr_dof )

        .field_readonly( "alpha_pt_mean", &bm::bvarw::alpha_pt_mean )
        .field_readonly( "alpha_pt_var", &bm::bvarw::alpha_pt_var )

        .field_readonly( "Sigma_pt_mean", &bm::bvarw::Sigma_pt_mean )
        .field_readonly( "Sigma_pt_dof", &bm::bvarw::Sigma_pt_dof )

        .field_readonly( "beta_draws", &bm::bvarw::beta_draws )
        .field_readonly( "Sigma_draws", &bm::bvarw::Sigma_draws )
    ;

    class_<bvarw_R>( "bvarw" )
        .derives<bm::bvarw>( "bvarw_cpp" )
        .default_constructor()

        .method( "build", build_1 )
        .method( "build", build_2 )
        .method( "reset_draws", &bvarw_R::reset_draws_R )
        .method( "prior", prior_3 )
        .method( "gibbs", &bvarw_R::gibbs_R )
        .method( "IRF", &bvarw_R::IRF_R )
        .method( "forecast", forecast_1 )
        .method( "forecast", forecast_2 )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

void bvarw_R::build_R(const arma::mat& data_raw, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarw_R::build_R(const arma::mat& data_raw, const arma::mat& data_ext, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,data_ext,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarw_R::reset_draws_R()
{
    try {
        this->reset_draws();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarw_R::prior_R(const arma::vec& coef_prior, double Xi_beta, double Xi_Sigma, int gamma)
{
    try {
        this->prior(coef_prior,Xi_beta,Xi_Sigma,gamma);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarw_R::gibbs_R(int n_draws, int n_burnin)
{
    try {
        this->gibbs(n_draws,n_burnin);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

SEXP bvarw_R::IRF_R(int n_irf_periods)
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

SEXP bvarw_R::forecast_R(int n_horizon, bool incl_shocks)
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

SEXP bvarw_R::forecast_R(const arma::mat& Y_T, int n_horizon, bool incl_shocks)
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
