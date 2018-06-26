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

RCPP_EXPOSED_CLASS(uhlig_R)
RCPP_EXPOSED_CLASS(dsge_uhlig_R)

RCPP_MODULE(dsgevar_uhlig_module)
{
    using namespace Rcpp ;
  
    class_<bm::dsgevar<bm::uhlig>>( "dsgevar_uhlig_cpp" )
        .default_constructor()

        // basic objects

        .field( "cons_term", &bm::dsgevar<bm::uhlig>::cons_term )
        .field( "p", &bm::dsgevar<bm::uhlig>::p )

        .field_readonly( "lambda", &bm::dsgevar<bm::uhlig>::lambda )
        .field_readonly( "c_int", &bm::dsgevar<bm::uhlig>::c_int )
        .field_readonly( "n", &bm::dsgevar<bm::uhlig>::n )
        .field_readonly( "M", &bm::dsgevar<bm::uhlig>::M )
        .field_readonly( "K", &bm::dsgevar<bm::uhlig>::K )

        .field_readonly( "Y", &bm::dsgevar<bm::uhlig>::Y )
        .field_readonly( "X", &bm::dsgevar<bm::uhlig>::X )

        .field_readonly( "alpha_pt_mean", &bm::dsgevar<bm::uhlig>::alpha_pt_mean )
        .field_readonly( "alpha_pt_var", &bm::dsgevar<bm::uhlig>::alpha_pt_var )
        .field_readonly( "Sigma_pt_mean", &bm::dsgevar<bm::uhlig>::Sigma_pt_mean )

        .field_readonly( "beta_draws", &bm::dsgevar<bm::uhlig>::beta_draws )
        .field_readonly( "Sigma_draws", &bm::dsgevar<bm::uhlig>::Sigma_draws )
    ;

    class_<dsgevar_uhlig_R>( "dsgevar_uhlig" )
        .derives<bm::dsgevar<bm::uhlig>>( "dsgevar_uhlig_cpp" )
        .default_constructor()

        .field( "opt_initial_lb", &dsgevar_uhlig_R::opt_initial_lb )
        .field( "opt_initial_ub", &dsgevar_uhlig_R::opt_initial_ub )
        .field( "mcmc_initial_lb", &dsgevar_uhlig_R::mcmc_initial_lb )
        .field( "mcmc_initial_ub", &dsgevar_uhlig_R::mcmc_initial_ub )

        .property( "dsge_draws", &dsgevar_uhlig_R::get_dsge_draws )
        .property( "lrem", &dsgevar_uhlig_R::get_lrem_R, &dsgevar_uhlig_R::set_lrem_R )
        .property( "dsge", &dsgevar_uhlig_R::get_dsge_R, &dsgevar_uhlig_R::set_dsge_R )

        .method( "build", &dsgevar_uhlig_R::build_R )

        .method( "get_model_fn", &dsgevar_uhlig_R::get_model_fn )
        .method( "set_model_fn", &dsgevar_uhlig_R::set_model_fn )
        .method( "eval_model", &dsgevar_uhlig_R::eval_model )

        .method( "set_bounds", &dsgevar_uhlig_R::set_bounds_R )
        .method( "set_prior", &dsgevar_uhlig_R::set_prior_R )

        .method( "estim_mode", &dsgevar_uhlig_R::estim_mode_R )
        .method( "estim_mcmc", &dsgevar_uhlig_R::estim_mcmc_R )

        .method( "mode_check", &dsgevar_uhlig_R::mode_check_R )

        .method( "IRF", &dsgevar_uhlig_R::IRF_R )
        .method( "forecast", &dsgevar_uhlig_R::forecast_R )
        .method( "states", &dsgevar_uhlig_R::state_filter_R )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

SEXP
dsgevar_uhlig_R::get_model_fn()
{
    return model_fn_SEXP;
}

void
dsgevar_uhlig_R::set_model_fn(SEXP model_fn_inp)
{
    try {
        if (TYPEOF(model_fn_inp) != CLOSXP) 
        {
            Rf_error("BMR: model_fn must be a function");
        }

        model_fn_SEXP = model_fn_inp;
        
        dsge_obj.model_fn = std::bind(&dsgevar_uhlig_R::model_fn_R,this,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::placeholders::_5,std::placeholders::_6);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

//

void
dsgevar_uhlig_R::eval_model(Rcpp::NumericVector pars_inp)
{
    try {

        Rcpp::Function pars_fn = Rcpp::as<Rcpp::Function>(model_fn_SEXP);
        Rcpp::List mats_out = pars_fn(pars_inp);

        //

        dsge_obj.lrem_obj.A = arma::mat( REAL(VECTOR_ELT(mats_out, 0)), Rf_nrows(VECTOR_ELT(mats_out, 0)), Rf_ncols(VECTOR_ELT(mats_out, 0)), false, true );
        dsge_obj.lrem_obj.B = arma::mat( REAL(VECTOR_ELT(mats_out, 1)), Rf_nrows(VECTOR_ELT(mats_out, 1)), Rf_ncols(VECTOR_ELT(mats_out, 1)), false, true );
        dsge_obj.lrem_obj.C = arma::mat( REAL(VECTOR_ELT(mats_out, 2)), Rf_nrows(VECTOR_ELT(mats_out, 2)), Rf_ncols(VECTOR_ELT(mats_out, 2)), false, true );
        dsge_obj.lrem_obj.D = arma::mat( REAL(VECTOR_ELT(mats_out, 3)), Rf_nrows(VECTOR_ELT(mats_out, 3)), Rf_ncols(VECTOR_ELT(mats_out, 3)), false, true );

        dsge_obj.lrem_obj.F = arma::mat( REAL(VECTOR_ELT(mats_out, 4)), Rf_nrows(VECTOR_ELT(mats_out, 4)), Rf_ncols(VECTOR_ELT(mats_out, 4)), false, true );
        dsge_obj.lrem_obj.G = arma::mat( REAL(VECTOR_ELT(mats_out, 5)), Rf_nrows(VECTOR_ELT(mats_out, 5)), Rf_ncols(VECTOR_ELT(mats_out, 5)), false, true );
        dsge_obj.lrem_obj.H = arma::mat( REAL(VECTOR_ELT(mats_out, 6)), Rf_nrows(VECTOR_ELT(mats_out, 6)), Rf_ncols(VECTOR_ELT(mats_out, 6)), false, true );

        dsge_obj.lrem_obj.J = arma::mat( REAL(VECTOR_ELT(mats_out, 7)), Rf_nrows(VECTOR_ELT(mats_out, 7)), Rf_ncols(VECTOR_ELT(mats_out, 7)), false, true );
        dsge_obj.lrem_obj.K = arma::mat( REAL(VECTOR_ELT(mats_out, 8)), Rf_nrows(VECTOR_ELT(mats_out, 8)), Rf_ncols(VECTOR_ELT(mats_out, 8)), false, true );
        dsge_obj.lrem_obj.L = arma::mat( REAL(VECTOR_ELT(mats_out, 9)), Rf_nrows(VECTOR_ELT(mats_out, 9)), Rf_ncols(VECTOR_ELT(mats_out, 9)), false, true );
        dsge_obj.lrem_obj.M = arma::mat( REAL(VECTOR_ELT(mats_out, 10)), Rf_nrows(VECTOR_ELT(mats_out, 10)), Rf_ncols(VECTOR_ELT(mats_out, 10)), false, true );
        dsge_obj.lrem_obj.N = arma::mat( REAL(VECTOR_ELT(mats_out, 11)), Rf_nrows(VECTOR_ELT(mats_out, 11)), Rf_ncols(VECTOR_ELT(mats_out, 11)), false, true );

    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void
dsgevar_uhlig_R::model_fn_R(const arma::vec& pars_inp, bm::uhlig& lrem_obj_inp, arma::mat& shocks_cov_out, arma::mat& C_out, arma::mat& H_out, arma::mat& R_out)
{
    try {

        Rcpp::NumericVector inp_vec(pars_inp.begin(),pars_inp.end());

        Rcpp::Function pars_fn = Rcpp::as<Rcpp::Function>(model_fn_SEXP);
        Rcpp::List mats_out = pars_fn(inp_vec);

        //

        lrem_obj_inp.A = arma::mat( REAL(VECTOR_ELT(mats_out, 0)), Rf_nrows(VECTOR_ELT(mats_out, 0)), Rf_ncols(VECTOR_ELT(mats_out, 0)), false, true );
        lrem_obj_inp.B = arma::mat( REAL(VECTOR_ELT(mats_out, 1)), Rf_nrows(VECTOR_ELT(mats_out, 1)), Rf_ncols(VECTOR_ELT(mats_out, 1)), false, true );
        lrem_obj_inp.C = arma::mat( REAL(VECTOR_ELT(mats_out, 2)), Rf_nrows(VECTOR_ELT(mats_out, 2)), Rf_ncols(VECTOR_ELT(mats_out, 2)), false, true );
        lrem_obj_inp.D = arma::mat( REAL(VECTOR_ELT(mats_out, 3)), Rf_nrows(VECTOR_ELT(mats_out, 3)), Rf_ncols(VECTOR_ELT(mats_out, 3)), false, true );

        lrem_obj_inp.F = arma::mat( REAL(VECTOR_ELT(mats_out, 4)), Rf_nrows(VECTOR_ELT(mats_out, 4)), Rf_ncols(VECTOR_ELT(mats_out, 4)), false, true );
        lrem_obj_inp.G = arma::mat( REAL(VECTOR_ELT(mats_out, 5)), Rf_nrows(VECTOR_ELT(mats_out, 5)), Rf_ncols(VECTOR_ELT(mats_out, 5)), false, true );
        lrem_obj_inp.H = arma::mat( REAL(VECTOR_ELT(mats_out, 6)), Rf_nrows(VECTOR_ELT(mats_out, 6)), Rf_ncols(VECTOR_ELT(mats_out, 6)), false, true );

        lrem_obj_inp.J = arma::mat( REAL(VECTOR_ELT(mats_out, 7)), Rf_nrows(VECTOR_ELT(mats_out, 7)), Rf_ncols(VECTOR_ELT(mats_out, 7)), false, true );
        lrem_obj_inp.K = arma::mat( REAL(VECTOR_ELT(mats_out, 8)), Rf_nrows(VECTOR_ELT(mats_out, 8)), Rf_ncols(VECTOR_ELT(mats_out, 8)), false, true );
        lrem_obj_inp.L = arma::mat( REAL(VECTOR_ELT(mats_out, 9)), Rf_nrows(VECTOR_ELT(mats_out, 9)), Rf_ncols(VECTOR_ELT(mats_out, 9)), false, true );
        lrem_obj_inp.M = arma::mat( REAL(VECTOR_ELT(mats_out, 10)), Rf_nrows(VECTOR_ELT(mats_out, 10)), Rf_ncols(VECTOR_ELT(mats_out, 10)), false, true );
        lrem_obj_inp.N = arma::mat( REAL(VECTOR_ELT(mats_out, 11)), Rf_nrows(VECTOR_ELT(mats_out, 11)), Rf_ncols(VECTOR_ELT(mats_out, 11)), false, true );

        shocks_cov_out = arma::mat( REAL(VECTOR_ELT(mats_out, 12)), Rf_nrows(VECTOR_ELT(mats_out, 12)), Rf_ncols(VECTOR_ELT(mats_out, 12)), false, true );

        C_out = arma::mat( REAL(VECTOR_ELT(mats_out, 13)), Rf_nrows(VECTOR_ELT(mats_out, 13)), Rf_ncols(VECTOR_ELT(mats_out, 13)), false, true );
        H_out = arma::mat( REAL(VECTOR_ELT(mats_out, 14)), Rf_nrows(VECTOR_ELT(mats_out, 14)), Rf_ncols(VECTOR_ELT(mats_out, 14)), false, true );
        R_out = arma::mat( REAL(VECTOR_ELT(mats_out, 15)), Rf_nrows(VECTOR_ELT(mats_out, 15)), Rf_ncols(VECTOR_ELT(mats_out, 15)), false, true );

    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

// get and set LREM

uhlig_R
dsgevar_uhlig_R::get_lrem_R()
{
    uhlig_R lrem_obj_out = static_cast<uhlig_R&>(dsge_obj.lrem_obj);

    return lrem_obj_out;
}

void
dsgevar_uhlig_R::set_lrem_R(uhlig_R lrem_obj_inp)
{
    try {
        dsge_obj.lrem_obj = static_cast<bm::uhlig&>(lrem_obj_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

// get and set DSGE model

dsge_uhlig_R
dsgevar_uhlig_R::get_dsge_R()
{
    dsge_uhlig_R dsge_obj_out = static_cast<dsge_uhlig_R&>(dsge_obj);

    return dsge_obj_out;
}

void
dsgevar_uhlig_R::set_dsge_R(dsge_uhlig_R dsge_obj_inp)
{
    try {
        dsge_obj = static_cast<bm::dsge<bm::uhlig>&>(dsge_obj_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

// build and access

void
dsgevar_uhlig_R::build_R(const arma::mat& data_raw, bool cons_term_inp, int p_inp, double lambda_inp)
{
    try {
        this->build(data_raw,cons_term_inp,p_inp,lambda_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

arma::mat
dsgevar_uhlig_R::get_dsge_draws()
{
    return dsge_obj.dsge_draws;
}

// set bounds and prior

void
dsgevar_uhlig_R::set_bounds_R(arma::vec lower_bounds_inp, arma::vec upper_bounds_inp)
{
    try {
        const int n_vals = lower_bounds_inp.n_elem;

        for (int i=0; i < n_vals; i++) {
            if (!R_FINITE(lower_bounds_inp[i])) { // R_NegInf
                lower_bounds_inp[i] = - bm::inf;
            }

            if (!R_FINITE(upper_bounds_inp[i])) { // R_PosInf
                upper_bounds_inp[i] = bm::inf;
            }
        }

        this->set_bounds(lower_bounds_inp,upper_bounds_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void
dsgevar_uhlig_R::set_prior_R(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp)
{
    try {
        this->set_prior(prior_form_inp,prior_pars_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

//

SEXP
dsgevar_uhlig_R::estim_mode_R(const arma::vec& initial_vals, bool calc_vcov)
{
    try {
        optim::algo_settings settings;

        if (opt_initial_lb.n_elem > 0) {
            settings.de_initial_lb = opt_initial_lb;
        }
        if (opt_initial_ub.n_elem > 0) {
            settings.de_initial_ub = opt_initial_ub;
        }

        settings.de_check_freq = 50;

        arma::vec res;
        arma::mat vcov_mat;

        if (calc_vcov) {
            res = this->estim_mode(initial_vals,&vcov_mat,&settings);
        } else {
            res = this->estim_mode(initial_vals,nullptr,&settings);
        }

        return Rcpp::List::create(Rcpp::Named("mode_vals") = res, Rcpp::Named("vcov_mat") = vcov_mat);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

void
dsgevar_uhlig_R::estim_mcmc_R(const arma::vec& initial_vals, int n_pop, int n_gen, int n_burnin)
{
    try {
        mcmc::algo_settings_t settings;

        if (mcmc_initial_lb.n_elem > 0) {
            settings.de_initial_lb = mcmc_initial_lb;
        }
        if (mcmc_initial_ub.n_elem > 0) {
            settings.de_initial_ub = mcmc_initial_ub;
        }

        settings.de_n_pop = n_pop;
        settings.de_n_gen = n_gen;
        settings.de_n_burnin = n_burnin;

        this->estim_mcmc(initial_vals,&settings);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

//

SEXP
dsgevar_uhlig_R::mode_check_R(const arma::vec& mode_vals, int grid_size, double scale_val)
{
    try {
        
        arma::cube mode_check_vals = bm::mode_check(*this,mode_vals,grid_size);

        return Rcpp::List::create(Rcpp::Named("mode_check_vals") = mode_check_vals);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

//

SEXP
dsgevar_uhlig_R::IRF_R(int n_irf_periods)
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

SEXP
dsgevar_uhlig_R::forecast_R(int n_horizon, bool incl_shocks)
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

SEXP
dsgevar_uhlig_R::state_filter_R()
{
    try {
        arma::cube states = this->state_filter();

        return Rcpp::List::create(Rcpp::Named("filter_vals") = states);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}
