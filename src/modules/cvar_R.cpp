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

#include "headers/var/cvar.hpp"

RCPP_MODULE(cvar_module)
{
    using namespace Rcpp ;

    // function overloading requires some trickery
    void (cvar::*data_1)(arma::mat) = &cvar::data ;
    void (cvar::*data_2)(arma::mat,arma::mat) = &cvar::data ;

    arma::cube (cvar::*forecast_1)(int, bool) = &cvar::forecast ;
    arma::cube (cvar::*forecast_2)(arma::mat, int, bool) = &cvar::forecast ;
  
    // now we can declare the class
    class_<cvar>( "R_cvar" )

    .default_constructor()

    // basic objects
    .field( "cons_term", &cvar::cons_term )
    .field( "p", &cvar::p )

    // read only objects
    .field_readonly( "n", &cvar::n )
    .field_readonly( "M", &cvar::M )
    .field_readonly( "K", &cvar::K )
    .field_readonly( "n_ext_vars", &cvar::n_ext_vars )

    .field_readonly( "beta_hat", &cvar::beta_hat )
    .field_readonly( "Sigma_hat", &cvar::Sigma_hat )

    .field_readonly( "beta_draws", &cvar::beta_draws )
    .field_readonly( "Sigma_draws", &cvar::Sigma_draws )
    .field_readonly( "irfs", &cvar::irfs )

    // member functions
    .method( "data", data_1 )
    .method( "data", data_2 )
    .method( "setup", &cvar::setup )
    .method( "boot", &cvar::boot )
    .method( "IRF", &cvar::IRF )
    .method( "forecast", forecast_1 )
    .method( "forecast", forecast_2 )
    ;
}
