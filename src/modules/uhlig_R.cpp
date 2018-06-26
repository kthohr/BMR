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

RCPP_MODULE(uhlig_module)
{
    using namespace Rcpp ;
  
    class_<bm::uhlig>( "uhlig_cpp" )
        .default_constructor()

        // basic objects
        .field( "A", &bm::uhlig::A )
        .field( "B", &bm::uhlig::B )
        .field( "C", &bm::uhlig::C )
        .field( "D", &bm::uhlig::D )

        .field( "F", &bm::uhlig::F )
        .field( "G", &bm::uhlig::G )
        .field( "H", &bm::uhlig::H )

        .field( "J", &bm::uhlig::J )
        .field( "K", &bm::uhlig::K )
        .field( "L", &bm::uhlig::L )
        .field( "M", &bm::uhlig::M )
        .field( "N", &bm::uhlig::N )

        .field( "shocks_cov", &bm::uhlig::shocks_cov )

        // read only objects
        .field_readonly( "P_sol", &bm::uhlig::P_sol )
        .field_readonly( "Q_sol", &bm::uhlig::Q_sol )
        .field_readonly( "R_sol", &bm::uhlig::R_sol )
        .field_readonly( "S_sol", &bm::uhlig::S_sol )

        .field_readonly( "F_state", &bm::uhlig::F_state )
        .field_readonly( "G_state", &bm::uhlig::G_state )
    ;

    class_<uhlig_R>( "uhlig" )
        .derives<bm::uhlig>( "uhlig_cpp" )
        .default_constructor()

        .method( "build", &uhlig_R::build_R )
        .method( "build_pre", &uhlig_R::build_pre_R )
        .method( "build_exp", &uhlig_R::build_exp_R )
        .method( "build_exog", &uhlig_R::build_exog_R )

        .method( "solve", &uhlig_R::solve_R )
        .method( "state_space", &uhlig_R::state_space_R )
        .method( "simulate", &uhlig_R::simulate_R )
        .method( "IRF", &uhlig_R::IRF_R )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

void uhlig_R::build_R(const arma::mat& A_inp, const arma::mat& B_inp, const arma::mat& C_inp, const arma::mat& D_inp, 
                      const arma::mat& F_inp, const arma::mat& G_inp, const arma::mat& H_inp, const arma::mat& J_inp, 
                      const arma::mat& K_inp, const arma::mat& L_inp, const arma::mat& M_inp, const arma::mat& N_inp)
{
    try {
        this->build(A_inp,B_inp,C_inp,D_inp,F_inp,G_inp,H_inp,J_inp,K_inp,L_inp,M_inp,N_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void uhlig_R::build_pre_R(const arma::mat& A_inp, const arma::mat& B_inp, const arma::mat& C_inp, const arma::mat& D_inp)
{
    try {
        this->build_pre(A_inp,B_inp,C_inp,D_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void uhlig_R::build_exp_R(const arma::mat& F_inp, const arma::mat& G_inp, const arma::mat& H_inp, const arma::mat& J_inp, 
                          const arma::mat& K_inp, const arma::mat& L_inp, const arma::mat& M_inp)
{
    try {
        this->build_exp(F_inp,G_inp,H_inp,J_inp,K_inp,L_inp,M_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void uhlig_R::build_exog_R(const arma::mat& N_inp)
{
    try {
        this->build_exog(N_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void uhlig_R::solve_R()
{
    try {
        this->solve();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void uhlig_R::state_space_R()
{
    try {
        this->state_space();
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

SEXP uhlig_R::simulate_R(int n_sim_periods, int n_burnin)
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

SEXP uhlig_R::IRF_R(int n_irf_periods)
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
