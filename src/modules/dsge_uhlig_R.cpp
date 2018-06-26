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

RCPP_MODULE(dsge_uhlig_module)
{
    using namespace Rcpp ;
  
    class_<bm::dsge<bm::uhlig>>( "dsge_uhlig_cpp" )
        .default_constructor()

        // basic objects
        .field( "estim_data", &bm::dsge<bm::uhlig>::estim_data )
        .field( "filter_choice", &bm::dsge<bm::uhlig>::filter_choice )

        .field_readonly( "prior_form", &bm::dsge<bm::uhlig>::prior_form )
        .field_readonly( "prior_pars", &bm::dsge<bm::uhlig>::prior_pars )

        .field_readonly( "lower_bounds", &bm::dsge<bm::uhlig>::lower_bounds )
        .field_readonly( "upper_bounds", &bm::dsge<bm::uhlig>::upper_bounds )

        .field_readonly( "dsge_draws", &bm::dsge<bm::uhlig>::dsge_draws )
    ;

    class_<dsge_uhlig_R>( "dsge_uhlig" )
        .derives<bm::dsge<bm::uhlig>>( "dsge_uhlig_cpp" )
        // .constructor<Rcpp::Function>()
        .default_constructor()

        .field( "opt_initial_lb", &dsge_uhlig_R::opt_initial_lb )
        .field( "opt_initial_ub", &dsge_uhlig_R::opt_initial_ub )
        .field( "mcmc_initial_lb", &dsge_uhlig_R::mcmc_initial_lb )
        .field( "mcmc_initial_ub", &dsge_uhlig_R::mcmc_initial_ub )

        .property( "lrem", &dsge_uhlig_R::get_lrem_R, &dsge_uhlig_R::set_lrem_R )
        // .property( "model", &dsge_uhlig_R::get_model_fn, &dsge_uhlig_R::set_model_fn )

        .method( "get_model_fn", &dsge_uhlig_R::get_model_fn )
        .method( "set_model_fn", &dsge_uhlig_R::set_model_fn )
        .method( "eval_model", &dsge_uhlig_R::eval_model )

        .method( "set_bounds", &dsge_uhlig_R::set_bounds_R )
        .method( "set_prior", &dsge_uhlig_R::set_prior_R )

        .method( "estim_mode", &dsge_uhlig_R::estim_mode_R )
        .method( "estim_mcmc", &dsge_uhlig_R::estim_mcmc_R )

        .method( "mode_check", &dsge_uhlig_R::mode_check_R )

        .method( "IRF", &dsge_uhlig_R::IRF_R )
        .method( "forecast", &dsge_uhlig_R::forecast_R )
        .method( "states", &dsge_uhlig_R::state_filter_R )
    ;
}

//
// wrapper functions to catch errors and handle memory pointers

SEXP
dsge_uhlig_R::get_model_fn()
{
    return model_fn_SEXP;
}

void
dsge_uhlig_R::set_model_fn(SEXP model_fn_inp)
{
    try {
        if (TYPEOF(model_fn_inp) != CLOSXP) 
        {
            Rf_error("BMR: model_fn must be a function");
        }

        model_fn_SEXP = model_fn_inp;

        // Rcpp::Function pars_fn = Rcpp::as<Rcpp::Function>(model_fn_SEXP);
        // model_fn_Rcpp = &pars_fn;
        // model_fn_Rcpp = reinterpret_cast<void*>(&pars_fn);

        // *(Rcpp::Function *)model_fn_Rcpp = pars_fn; 
        
        model_fn = std::bind(&dsge_uhlig_R::model_fn_R,this,std::placeholders::_1,std::placeholders::_2,std::placeholders::_3,std::placeholders::_4,std::placeholders::_5,std::placeholders::_6);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

//

void
dsge_uhlig_R::eval_model(Rcpp::NumericVector pars_inp)
{
    try {

        Rcpp::Function pars_fn = Rcpp::as<Rcpp::Function>(model_fn_SEXP);
        Rcpp::List mats_out = pars_fn(pars_inp);

        //

        lrem_obj.A = arma::mat( REAL(VECTOR_ELT(mats_out, 0)), Rf_nrows(VECTOR_ELT(mats_out, 0)), Rf_ncols(VECTOR_ELT(mats_out, 0)), false, true );
        lrem_obj.B = arma::mat( REAL(VECTOR_ELT(mats_out, 1)), Rf_nrows(VECTOR_ELT(mats_out, 1)), Rf_ncols(VECTOR_ELT(mats_out, 1)), false, true );
        lrem_obj.C = arma::mat( REAL(VECTOR_ELT(mats_out, 2)), Rf_nrows(VECTOR_ELT(mats_out, 2)), Rf_ncols(VECTOR_ELT(mats_out, 2)), false, true );
        lrem_obj.D = arma::mat( REAL(VECTOR_ELT(mats_out, 3)), Rf_nrows(VECTOR_ELT(mats_out, 3)), Rf_ncols(VECTOR_ELT(mats_out, 3)), false, true );

        lrem_obj.F = arma::mat( REAL(VECTOR_ELT(mats_out, 4)), Rf_nrows(VECTOR_ELT(mats_out, 4)), Rf_ncols(VECTOR_ELT(mats_out, 4)), false, true );
        lrem_obj.G = arma::mat( REAL(VECTOR_ELT(mats_out, 5)), Rf_nrows(VECTOR_ELT(mats_out, 5)), Rf_ncols(VECTOR_ELT(mats_out, 5)), false, true );
        lrem_obj.H = arma::mat( REAL(VECTOR_ELT(mats_out, 6)), Rf_nrows(VECTOR_ELT(mats_out, 6)), Rf_ncols(VECTOR_ELT(mats_out, 6)), false, true );

        lrem_obj.J = arma::mat( REAL(VECTOR_ELT(mats_out, 7)), Rf_nrows(VECTOR_ELT(mats_out, 7)), Rf_ncols(VECTOR_ELT(mats_out, 7)), false, true );
        lrem_obj.K = arma::mat( REAL(VECTOR_ELT(mats_out, 8)), Rf_nrows(VECTOR_ELT(mats_out, 8)), Rf_ncols(VECTOR_ELT(mats_out, 8)), false, true );
        lrem_obj.L = arma::mat( REAL(VECTOR_ELT(mats_out, 9)), Rf_nrows(VECTOR_ELT(mats_out, 9)), Rf_ncols(VECTOR_ELT(mats_out, 9)), false, true );
        lrem_obj.M = arma::mat( REAL(VECTOR_ELT(mats_out, 10)), Rf_nrows(VECTOR_ELT(mats_out, 10)), Rf_ncols(VECTOR_ELT(mats_out, 10)), false, true );
        lrem_obj.N = arma::mat( REAL(VECTOR_ELT(mats_out, 11)), Rf_nrows(VECTOR_ELT(mats_out, 11)), Rf_ncols(VECTOR_ELT(mats_out, 11)), false, true );

    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

void
dsge_uhlig_R::model_fn_R(const arma::vec& pars_inp, bm::uhlig& lrem_obj_inp, arma::mat& shocks_cov_out, arma::mat& C_out, arma::mat& H_out, arma::mat& R_out)
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
dsge_uhlig_R::get_lrem_R()
{
    uhlig_R lrem_obj_out = static_cast<uhlig_R&>(lrem_obj);

    return lrem_obj_out;
}

void
dsge_uhlig_R::set_lrem_R(uhlig_R lrem_obj_inp)
{
    try {
        lrem_obj = static_cast<bm::uhlig&>(lrem_obj_inp);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
}

// set bounds and prior

void
dsge_uhlig_R::set_bounds_R(arma::vec lower_bounds_inp, arma::vec upper_bounds_inp)
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
dsge_uhlig_R::set_prior_R(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp)
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
dsge_uhlig_R::estim_mode_R(const arma::vec& initial_vals, bool calc_vcov)
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
dsge_uhlig_R::estim_mcmc_R(const arma::vec& initial_vals, int n_pop, int n_gen, int n_burnin)
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
dsge_uhlig_R::mode_check_R(const arma::vec& mode_vals, int grid_size, double scale_val)
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
dsge_uhlig_R::IRF_R(int n_irf_periods, bool observ_irfs)
{
    try {
        arma::cube irf_vals = this->IRF(n_irf_periods,observ_irfs);
        return Rcpp::List::create(Rcpp::Named("irf_vals") = irf_vals);
    } catch( std::exception &ex ) {
        forward_exception_to_r( ex );
    } catch(...) {
        ::Rf_error( "BMR: C++ exception (unknown reason)" );
    }
    return R_NilValue;
}

SEXP
dsge_uhlig_R::forecast_R(int n_horizon, bool incl_shocks)
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
dsge_uhlig_R::state_filter_R()
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
