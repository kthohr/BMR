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

RCPP_MODULE(bvarm_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (bvarm_R::*build_1)(const arma::mat&, bool, int) = &bvarm_R::build_R;
    void (bvarm_R::*build_2)(const arma::mat&,const arma::mat&, bool, int) = &bvarm_R::build_R;

    void (bvarm_R::*prior_1)(const arma::vec&) = &bvarm_R::prior_R;
    void (bvarm_R::*prior_2)(const arma::vec&, int, int, double, double, double, double) = &bvarm_R::prior_R;

    SEXP (bvarm_R::*forecast_1)(int, bool) = &bvarm_R::forecast_R;
    SEXP (bvarm_R::*forecast_2)(const arma::mat&, int, bool) = &bvarm_R::forecast_R;
  
    // now we can declare the class
    class_<bm::bvarm>( "bvarm_cpp" )
        .default_constructor()

        // basic objects
        .field( "cons_term", &bm::bvarm::cons_term )
        .field( "p", &bm::bvarm::p )

        .field( "only_stationary_draws", &bm::bvarm::only_stationary_draws )
        .field( "irfs_lr_restrict", &bm::bvarm::irfs_lr_restrict )

        // read only objects
        .field_readonly( "var_type", &bm::bvarm::var_type )
        .field_readonly( "decay_type", &bm::bvarm::decay_type )

        .field_readonly( "c_int", &bm::bvarm::c_int )
        .field_readonly( "n", &bm::bvarm::n )
        .field_readonly( "M", &bm::bvarm::M )
        .field_readonly( "K", &bm::bvarm::K )
        .field_readonly( "n_ext_vars", &bm::bvarm::n_ext_vars )

        .field_readonly( "Y", &bm::bvarm::Y )
        .field_readonly( "X", &bm::bvarm::X )

        .field_readonly( "alpha_hat", &bm::bvarm::alpha_hat )
        .field_readonly( "Sigma_hat", &bm::bvarm::Sigma_hat )

        .field_readonly( "alpha_pr_mean", &bm::bvarm::alpha_pr_mean )
        .field_readonly( "alpha_pr_var", &bm::bvarm::alpha_pr_var )

        .field_readonly( "alpha_pt_mean", &bm::bvarm::alpha_pt_mean )
        .field_readonly( "alpha_pt_var", &bm::bvarm::alpha_pt_var )

        .field_readonly( "beta_draws", &bm::bvarm::beta_draws )
    ;

    class_<bvarm_R>( "bvarm" )
        .derives<bm::bvarm>( "bvarm_cpp" )
        .default_constructor()

        .method( "build", build_1 )
        .method( "build", build_2 )
        .method( "reset_draws", &bvarm_R::reset_draws_R )
        .method( "prior", prior_1 )
        .method( "prior", prior_2 )
        .method( "gibbs", &bvarm_R::gibbs_R )
        .method( "IRF", &bvarm_R::IRF_R )
        .method( "FEVD", &bvarm_R::FEVD_R )
        .method( "forecast", forecast_1 )
        .method( "forecast", forecast_2 )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

void bvarm_R::build_R(const arma::mat& data_raw, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarm_R::build_R(const arma::mat& data_raw, const arma::mat& data_ext, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,data_ext,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarm_R::reset_draws_R()
{
    try {
        this->reset_draws();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarm_R::prior_R(const arma::vec& coef_prior)
{
    try {
        this->prior(coef_prior,1,1,0.5,0.5,1.0,2.0);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarm_R::prior_R(const arma::vec& coef_prior, int var_type_inp, int decay_type_inp, double HP_1_inp, double HP_2_inp, double HP_3_inp, double HP_4_inp)
{
    try {
        this->prior(coef_prior,var_type_inp,decay_type_inp,HP_1_inp,HP_2_inp,HP_3_inp,HP_4_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvarm_R::gibbs_R(int n_draws)
{
    try {
        this->gibbs(n_draws);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

SEXP bvarm_R::IRF_R(int n_irf_periods)
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

SEXP bvarm_R::FEVD_R(int n_periods)
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

SEXP bvarm_R::forecast_R(int n_horizon, bool incl_shocks)
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

SEXP bvarm_R::forecast_R(const arma::mat& Y_T, int n_horizon, bool incl_shocks)
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
