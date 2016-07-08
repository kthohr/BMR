/*################################################################################
  ##
  ##   R package BMR by Keith O'Hara Copyright (C) 2011-2016
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
  ################################################################################*/

#include <RcppArmadillo.h>

#include "BMR_misc.hpp"
#include "BMR_stats.hpp"

#include "headers/bvar/bvarw.hpp"

RCPP_MODULE(bvarw_module){
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (bvarw::*data_1)(arma::mat) = &bvarw::data ;
    void (bvarw::*data_2)(arma::mat,arma::mat) = &bvarw::data ;

    void (bvarw::*prior_1)(arma::vec, double, arma::mat, int) = &bvarw::prior ;
    void (bvarw::*prior_2)(arma::vec, arma::mat, arma::mat, int) = &bvarw::prior ;

    arma::cube (bvarw::*forecast_1)(int, bool) = &bvarw::forecast ;
    arma::cube (bvarw::*forecast_2)(arma::mat, int, bool) = &bvarw::forecast ;
  
    // now we can declare the class
    class_<bvarw>( "R_bvarw" )

    .default_constructor()

    // basic objects
    .field( "cons_term", &bvarw::cons_term )
    .field( "p", &bvarw::p )

    // read only objects
    .field_readonly( "n", &bvarw::n )
    .field_readonly( "M", &bvarw::M )
    .field_readonly( "K", &bvarw::K )
    .field_readonly( "n_ext_vars", &bvarw::n_ext_vars )

    .field_readonly( "alpha_pr_mean", &bvarw::alpha_pr_mean )
    .field_readonly( "alpha_pr_var", &bvarw::alpha_pr_var )

    .field_readonly( "Sigma_pr_scale", &bvarw::Sigma_pr_scale )
    .field_readonly( "Sigma_pr_dof", &bvarw::Sigma_pr_dof )

    .field_readonly( "alpha_hat", &bvarw::alpha_hat )
    .field_readonly( "Sigma_hat", &bvarw::Sigma_hat )

    .field_readonly( "alpha_pt_mean", &bvarw::alpha_pt_mean )
    .field_readonly( "alpha_pt_var", &bvarw::alpha_pt_var )

    .field_readonly( "Sigma_pt_scale", &bvarw::Sigma_pt_scale )
    .field_readonly( "Sigma_pt_dof", &bvarw::Sigma_pt_dof )
    .field_readonly( "Sigma_pt_mean", &bvarw::Sigma_pt_mean )

    .field_readonly( "beta_draws", &bvarw::beta_draws )
    .field_readonly( "Sigma_draws", &bvarw::Sigma_draws )
    .field_readonly( "irfs", &bvarw::irfs )

    // member functions
    .method( "data", data_1 )
    .method( "data", data_2 )
    .method( "prior", prior_1 )
    .method( "prior", prior_2 )
    .method( "gibbs", &bvarw::gibbs )
    .method( "IRF", &bvarw::IRF )
    .method( "forecast", forecast_1 )
    .method( "forecast", forecast_2 )
    ;
}
