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

RCPP_MODULE(bvartvp_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery

    // SEXP (bvartvp_R::*forecast_1)(int, bool) = &bvartvp_R::forecast_R;
    // SEXP (bvartvp_R::*forecast_2)(const arma::mat&, int, bool) = &bvartvp_R::forecast_R;
  
    // now we can declare the class
    class_<bm::bvartvp>( "bvartvp_cpp" )
        .default_constructor()

        // basic objects
        .field( "cons_term", &bm::bvartvp::cons_term )
        .field( "p", &bm::bvartvp::p )

        // read only objects
        .field_readonly( "c_int", &bm::bvartvp::c_int )
        .field_readonly( "n", &bm::bvartvp::n )
        .field_readonly( "M", &bm::bvartvp::M )
        .field_readonly( "K", &bm::bvartvp::K )
        .field_readonly( "n_ext_vars", &bm::bvartvp::n_ext_vars )
        .field_readonly( "tau", &bm::bvartvp::tau )

        .field_readonly( "alpha_hat", &bm::bvartvp::alpha_hat )
        .field_readonly( "Sigma_hat", &bm::bvartvp::Sigma_hat )

        .field_readonly( "alpha_pr_mean", &bm::bvartvp::alpha_pr_mean )
        .field_readonly( "alpha_pr_var", &bm::bvartvp::alpha_pr_var )

        .field_readonly( "Sigma_pr_scale", &bm::bvartvp::Sigma_pr_scale )
        .field_readonly( "Sigma_pr_dof", &bm::bvartvp::Sigma_pr_dof )

        .field_readonly( "alpha_pt_mean", &bm::bvartvp::alpha_pt_mean )

        .field_readonly( "Sigma_pt_mean", &bm::bvartvp::Sigma_pt_mean )
        .field_readonly( "Sigma_pt_dof", &bm::bvartvp::Sigma_pt_dof )

        .field_readonly( "Q_pt_mean", &bm::bvartvp::Q_pt_mean )
        .field_readonly( "Q_pt_dof", &bm::bvartvp::Q_pt_dof )

        .field_readonly( "alpha_draws", &bm::bvartvp::alpha_draws )
        .field_readonly( "Q_draws", &bm::bvartvp::Q_draws )
        .field_readonly( "Sigma_draws", &bm::bvartvp::Sigma_draws )
    ;

    class_<bvartvp_R>( "bvartvp" )
        .derives<bm::bvartvp>( "bvartvp_cpp" )
        .default_constructor()

        .method( "build", &bvartvp_R::build_R )
        .method( "prior", &bvartvp_R::prior_R )
        .method( "reset_draws", &bvartvp_R::reset_draws_R )
        .method( "gibbs", &bvartvp_R::gibbs_R )
        .method( "IRF", &bvartvp_R::IRF_R )
        // .method( "forecast", forecast_1 )
        // .method( "forecast", forecast_2 )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

void bvartvp_R::build_R(const arma::mat& data_raw, bool cons_term_inp, int p_inp)
{
    try {
        this->build(data_raw,cons_term_inp,p_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvartvp_R::reset_draws_R()
{
    try {
        this->reset_draws();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvartvp_R::prior_R(int tau_inp, double Xi_beta, double Xi_Q_inp, int gamma_Q, double Xi_Sigma, int gamma_S)
{
    try {
        this->prior(tau_inp,Xi_beta,Xi_Q_inp,gamma_Q,Xi_Sigma,gamma_S);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void bvartvp_R::gibbs_R(int n_draws, int n_burnin)
{
    try {
        this->gibbs(n_draws,n_burnin);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

SEXP bvartvp_R::IRF_R(int n_irf_periods, int time_period)
{
    try {
        arma::cube irf_vals = this->IRF(n_irf_periods,time_period);

        return Rcpp::List::create(Rcpp::Named("irf_vals") = irf_vals);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
