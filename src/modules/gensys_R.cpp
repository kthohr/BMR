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
#include "dsge_R.hpp"

RCPP_MODULE(gensys_module)
{
    using namespace Rcpp ;
  
    class_<bm::gensys>( "gensys_cpp" )
        .default_constructor()

        // basic objects
        .field( "Gamma_0", &bm::gensys::Gamma_0 )
        .field( "Gamma_1", &bm::gensys::Gamma_1 )
        .field( "Gamma_C", &bm::gensys::Gamma_C )
        .field( "Psi", &bm::gensys::Psi )
        .field( "Pi", &bm::gensys::Pi )

        .field( "shocks_cov", &bm::gensys::shocks_cov )

        // read only objects
        .field_readonly( "G_sol", &bm::gensys::G_sol )
        .field_readonly( "impact_sol", &bm::gensys::impact_sol )
        .field_readonly( "cons_sol", &bm::gensys::cons_sol )

        .field_readonly( "F_state", &bm::gensys::F_state )
        .field_readonly( "G_state", &bm::gensys::G_state )
    ;

    class_<gensys_R>( "gensys" )
        .derives<bm::gensys>( "gensys_cpp" )
        .default_constructor()

        .method( "build", &gensys_R::build_R )
        .method( "solve", &gensys_R::solve_R )
        .method( "state_space", &gensys_R::state_space_R )
        .method( "simulate", &gensys_R::simulate_R )
        .method( "IRF", &gensys_R::IRF_R )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

void gensys_R::build_R(const arma::mat& Gamma_0_inp, const arma::mat& Gamma_1_inp, const arma::mat& Gamma_C_inp, const arma::mat& Psi_inp, const arma::mat& Pi_inp)
{
    try {
        this->build(Gamma_0_inp,Gamma_1_inp,Gamma_C_inp,Psi_inp,Pi_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void gensys_R::solve_R()
{
    try {
        this->solve();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void gensys_R::state_space_R()
{
    try {
        this->state_space();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

SEXP gensys_R::simulate_R(int n_sim_periods, int n_burnin)
{
    try {
        arma::mat ret_mat = this->simulate(n_sim_periods,n_burnin);
        return Rcpp::List::create(Rcpp::Named("sim_vals") = ret_mat);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP gensys_R::IRF_R(int n_irf_periods)
{
    try {
        arma::cube ret_cube = this->IRF(n_irf_periods);
        return Rcpp::List::create(Rcpp::Named("irf_vals") = ret_cube);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
