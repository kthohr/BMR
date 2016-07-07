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

#include "headers/misc/embed.hpp"
#include "headers/stats/rmvnorm.hpp"
#include "headers/bvar/bvarm.hpp"

RCPP_MODULE(bvarm_module){
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (bvarm::*data_1)(arma::mat) = &bvarm::data ;
    void (bvarm::*data_2)(arma::mat,arma::mat) = &bvarm::data ;

    void (bvarm::*IRF_1)(int) = &bvarm::IRF ;
    void (bvarm::*IRF_2)(arma::mat,int) = &bvarm::IRF ;

    arma::cube (bvarm::*forecast_1)(arma::mat, int, bool) = &bvarm::forecast ;
    arma::cube (bvarm::*forecast_2)(arma::mat, int, bool, arma::mat) = &bvarm::forecast ;
  
    // now we can declare the class
    class_<bvarm>( "R_bvarm" )

    .default_constructor()

    // basic objects
    .field( "cons_term", &bvarm::cons_term )
    .field( "p", &bvarm::p )

    .field( "var_type", &bvarm::var_type )
    .field( "decay_type", &bvarm::decay_type )
    .field( "hyper_pars", &bvarm::hyper_pars )

    // read only objects
    .field_readonly( "n", &bvarm::n )
    .field_readonly( "M", &bvarm::M )
    .field_readonly( "K", &bvarm::K )
    .field_readonly( "n_ext_vars", &bvarm::n_ext_vars )

    .field_readonly( "alpha_pr_mean", &bvarm::alpha_pr_mean )
    .field_readonly( "alpha_pr_var", &bvarm::alpha_pr_var )

    .field_readonly( "alpha_hat", &bvarm::alpha_hat )
    .field_readonly( "Sigma", &bvarm::Sigma )

    .field_readonly( "alpha_pt_mean", &bvarm::alpha_pt_mean )
    .field_readonly( "alpha_pt_var", &bvarm::alpha_pt_var )

    .field_readonly( "beta_draws", &bvarm::beta_draws )
    .field_readonly( "irfs", &bvarm::irfs )

    // member functions
    .method( "data", data_1 )
    .method( "data", data_2 )
    .method( "prior", &bvarm::prior )
    .method( "gibbs", &bvarm::gibbs )
    .method( "IRF", IRF_1 )
    .method( "IRF", IRF_2 )
    .method( "forecast", forecast_1 )
    .method( "forecast", forecast_2 )
    ;
}
