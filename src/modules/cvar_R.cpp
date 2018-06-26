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

RCPP_MODULE(cvar_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (cvar_R::*build_1)(const arma::mat&,bool,int) = &cvar_R::build_R;
    void (cvar_R::*build_2)(const arma::mat&,const arma::mat&,bool,int) = &cvar_R::build_R;

    SEXP (cvar_R::*forecast_1)(int, bool) = &cvar_R::forecast_R;
    SEXP (cvar_R::*forecast_2)(const arma::mat&, int, bool) = &cvar_R::forecast_R;
  
    // now we can declare the class
    class_<bm::cvar>( "cvar_cpp" )
        .default_constructor()

        // basic objects
        .field( "cons_term", &bm::cvar::cons_term )
        .field( "p", &bm::cvar::p )

        .field( "irfs_lr_restrict", &bm::cvar::irfs_lr_restrict )

        // read only objects
        .field_readonly( "c_int", &bm::cvar::c_int )
        .field_readonly( "n", &bm::cvar::n )
        .field_readonly( "M", &bm::cvar::M )
        .field_readonly( "K", &bm::cvar::K )
        .field_readonly( "n_ext_vars", &bm::cvar::n_ext_vars )

        .field_readonly( "Y", &bm::cvar::Y )
        .field_readonly( "X", &bm::cvar::X )

        .field_readonly( "beta_hat", &bm::cvar::beta_hat )
        .field_readonly( "Sigma_hat", &bm::cvar::Sigma_hat )

        .field_readonly( "beta_draws", &bm::cvar::beta_draws )
        .field_readonly( "Sigma_draws", &bm::cvar::Sigma_draws )
    ;

    class_<cvar_R>( "cvar" )
        .derives<bm::cvar>( "cvar_cpp" )
        .default_constructor()

        .method( "build", build_1 )
        .method( "build", build_2 )
        .method( "reset_draws", &cvar_R::reset_draws_R )
        .method( "estim", &cvar_R::estim_R )
        .method( "boot", &cvar_R::boot_R )
        .method( "IRF", &cvar_R::IRF_R )
        .method( "FEVD", &cvar_R::FEVD_R )
        .method( "forecast", forecast_1 )
        .method( "forecast", forecast_2 )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

void cvar_R::build_R(const arma::mat& data_raw, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void cvar_R::build_R(const arma::mat& data_raw, const arma::mat& data_ext, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,data_ext,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void cvar_R::reset_draws_R()
{
    try {
        this->reset_draws();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void cvar_R::estim_R()
{
    try {
        this->estim();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void cvar_R::boot_R(int n_draws)
{
    try {
        this->boot(n_draws);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

SEXP cvar_R::IRF_R(int n_irf_periods)
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

SEXP cvar_R::FEVD_R(int n_periods)
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

SEXP cvar_R::forecast_R(int n_horizon, bool incl_shocks)
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

SEXP cvar_R::forecast_R(const arma::mat& Y_T, int n_horizon, bool incl_shocks)
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
