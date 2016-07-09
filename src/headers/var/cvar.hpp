/*################################################################################
  ##
  ##   MacroLib by Keith O'Hara Copyright (C) 2011-2016
  ##   This file is part of the C++ MacroLib library.
  ##
  ##   MacroLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   MacroLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/*
 * CVAR
 */

//#include "embed.hpp"

// cvar class
class cvar
{
    public:
        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags
        int M;          // number of endogenous variables
        int K;          // number of coefficients in each 
        int n_ext_vars; // number of 'external' variables

        arma::mat Y; // Y = X beta + e
        arma::mat X;

        arma::mat beta_hat;       // OLS estimate of beta
        arma::mat Sigma_hat;      // OLS-based estimation of covariance matrix of 'e'

        arma::cube beta_draws;    // bootstrap draws of beta
        arma::cube Sigma_draws;   // bootstrap draws of beta
        arma::cube irfs;          // irfs based on the bootstrap draws

        // member functions
        void data(arma::mat data_raw);
        void data(arma::mat data_raw, arma::mat data_ext);
        void setup();
        void boot(int n_draws);
        void IRF(int n_irf_periods);
        arma::cube forecast(int horizon, bool incl_shocks);
        arma::cube forecast(arma::mat Y_T, int horizon, bool incl_shocks);

    //private:
};

// This is broken into two functions (instead of using data(arma::mat data_raw, arma::mat* data_ext))
// because Rcpp modules can't handle pointers
void cvar::data(arma::mat data_raw)
{
    // data_raw is a n x M matrix with 'endogenous' variables
    c_int = !!cons_term;

    n = data_raw.n_rows;
    M = data_raw.n_cols;

    n_ext_vars = 0;
    K = c_int + M*p + n_ext_vars;
    //
    arma::mat data_lagged = embed(data_raw,p);
    //
    Y = data_lagged.cols(0,M-1);

    if (cons_term) {  // setup lagged variables
        X = arma::join_rows(arma::ones(n-p),data_lagged.cols(M,data_lagged.n_cols-1));
    } else {
        X = data_lagged.cols(M,data_lagged.n_cols-1);
    }
    //
}

void cvar::data(arma::mat data_raw, arma::mat data_ext)
{
    // data_raw is a n x M matrix with 'endogenous' variables
    // data_ext is a n x n_ext_vars matrix with 'exogenous' variables
    c_int = !!cons_term;

    n = data_raw.n_rows;
    M = data_raw.n_cols;

    n_ext_vars = data_ext.n_cols;
    K = c_int + M*p + n_ext_vars;
    //
    arma::mat data_lagged = embed(data_raw,p);
    //
    Y = data_lagged.cols(0,M-1);

    if (cons_term) {  // setup lagged variables
        X = arma::join_rows(arma::ones(n-p),data_lagged.cols(M,data_lagged.n_cols-1));
    } else {
        X = data_lagged.cols(M,data_lagged.n_cols-1);
    }
    //
    if (n_ext_vars > 0) {
        arma::mat X_ext = data_ext;
        X_ext.shed_rows(0,p-1);

        X = arma::join_rows(X,X_ext);
    }
    //
}

void cvar::setup()
{
    //
    beta_hat = arma::solve(X.t()*X,X.t()*Y);
    
    arma::mat epsilon = Y - X*beta_hat;
    Sigma_hat = epsilon.t() * epsilon / ((double) n - p); // MLE (not bias-corrected)
    //
}

void cvar::boot(int n_draws)
{
    int i,j;
    int K_adj = K - n_ext_vars;

    beta_draws.set_size(K, M, n_draws);
    Sigma_draws.set_size(M, M, n_draws);
    
    arma::mat beta_b  = beta_hat;
    arma::mat Sigma_b = Sigma_hat;
    //
    arma::ivec sampling_vec = arma::randi(n - p, arma::distr_param(0,n-p-1));
    arma::uvec eps_sample = arma::conv_to<arma::uvec>::from(sampling_vec);

    arma::mat epsilon_hat = Y - X * beta_hat;
    arma::mat epsilon_b = epsilon_hat.rows(eps_sample);

    arma::mat Y_b = Y;
    arma::mat X_b = X;

    arma::mat X_t_block;
    //
    for (i=0; i < n_draws; i++) {
        sampling_vec = arma::randi(n - p, arma::distr_param(0,n-p-1));
        eps_sample = arma::conv_to<arma::uvec>::from(sampling_vec);
        epsilon_b = epsilon_hat.rows(eps_sample);
        
        X_t_block = X.row(0);
        //
        for (j = 0; j < n-p; j++) {
            X_b.row(j) = X_t_block;
            Y_b.row(j) = X_t_block * beta_hat + epsilon_b.row(j);
            //
            if (K_adj > M + c_int) {
                X_t_block(0,arma::span(M+c_int,K_adj-1)) = X_t_block(0,arma::span(c_int,K_adj-M-1));
            }

            if (n_ext_vars > 0) {
                X_t_block(0,arma::span(K_adj,K-1)) = X(j,arma::span(K_adj,K-1)); // update external variables
            }

            X_t_block(0,arma::span(c_int,M-1+c_int)) = Y_b.row(j);
        }
        //
        beta_b = arma::solve(X_b.t()*X_b,X_b.t()*Y_b);
        
        epsilon_b = Y_b - X_b * beta_b;
        Sigma_b = epsilon_b.t()*epsilon_b / ((double) n-p); // MLE (not bias-corrected)
        //
        beta_draws.slice(i)  = beta_b;
        Sigma_draws.slice(i) = Sigma_b;
    }
    //
}

void cvar::IRF(int n_irf_periods)
{
    int i,j;

    int n_draws = beta_draws.n_slices;
    int K_adj = K - n_ext_vars;

    arma::mat Sigma_b, impact_mat;
    
    irfs.set_size(M, M, n_irf_periods*n_draws);
    //
    arma::mat beta_b(K_adj-c_int,M);        // b'th draw, minus coefficients on any external variables 
    arma::mat beta_b_trans(K_adj-c_int,M); 

    arma::mat impact_mat_b(K_adj-c_int,M);
    arma::mat impact_mat_h(M,M);
    //
    for (j=1; j<=n_draws; j++) {
        beta_b = beta_draws(arma::span(c_int,K_adj-1),arma::span(),arma::span(j-1,j-1));
        beta_b_trans = beta_b.t();

        Sigma_b = Sigma_draws.slice(j-1);
        impact_mat = arma::trans(arma::chol(Sigma_b));
        //
        impact_mat_b.zeros();
        impact_mat_b.rows(0,M-1) = impact_mat;

        irfs.slice((j-1)*n_irf_periods) = impact_mat;

        for (i=2; i<=n_irf_periods; i++) {
            impact_mat_h = beta_b_trans*impact_mat_b;
            irfs.slice((j-1)*n_irf_periods + (i-1)) = impact_mat_h;

            if(K_adj > M+c_int){
                impact_mat_b.rows(M,K_adj-c_int-1) = impact_mat_b.rows(0,K_adj-M-c_int-1);
            }

            impact_mat_b.rows(0,M-1) = impact_mat_h;
        }
    }
    //
}

arma::cube cvar::forecast(int horizon, bool incl_shocks)
{
    int i,j;

    int n_draws = beta_draws.n_slices;
    int K_adj = K - n_ext_vars;
    
    arma::mat beta_b(K_adj,M), Sigma_b(M,M);       // bth draw
    
    arma::mat Y_T  = arma::join_rows(Y.row(Y.n_rows-1),X(X.n_rows-1,arma::span(0,K_adj-M-1)));
    arma::mat Y_Th = Y_T;

    arma::mat Y_forecast(horizon,M);
    arma::cube forecast_mat(horizon, M, n_draws);
    //
    if (incl_shocks) {
        arma::mat chol_shock_cov = arma::trans(arma::chol(Sigma_hat));

        for (i=0; i<n_draws; i++) {
            beta_b  = beta_draws.slice(i);
            Sigma_b = Sigma_draws.slice(i);

            chol_shock_cov = arma::trans(arma::chol(Sigma_b));

            Y_forecast.zeros();
            Y_Th = Y_T;
            
            for (j=1; j<=horizon; j++) {
                Y_forecast.row(j-1) = Y_Th*beta_b + arma::trans(rmvnorm(chol_shock_cov,true));

                if (K_adj > M + c_int) {
                    Y_Th(0,arma::span(M+c_int,K_adj-1)) = Y_Th(0,arma::span(c_int,K_adj-M-1));
                }

                Y_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
            }
            //
            forecast_mat.slice(i) = Y_forecast;
        }
    } else {
        for (i=0; i<n_draws; i++) {
            beta_b = beta_draws.slice(i);

            Y_forecast.zeros();
            Y_Th = Y_T;

            for (j=1; j<=horizon; j++) {
                Y_forecast.row(j-1) = Y_Th*beta_b;

                if (K_adj > M + c_int) {
                    Y_Th(0,arma::span(M+c_int,K_adj-1)) = Y_Th(0,arma::span(c_int,K_adj-M-1));
                }

                Y_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
            }
            //
            forecast_mat.slice(i) = Y_forecast;
        }
    }
    //
    return forecast_mat;
}

arma::cube cvar::forecast(arma::mat Y_T, int horizon, bool incl_shocks)
{
    int i,j;

    int n_draws = beta_draws.n_slices;
    int K_adj = K - n_ext_vars;
    
    arma::mat beta_b(K_adj,M), Sigma_b(M,M);       // bth draw
    
    arma::mat Y_Th = Y_T;

    arma::mat Y_forecast(horizon,M);
    arma::cube forecast_mat(horizon, M, n_draws);
    //
    if (incl_shocks) {
        arma::mat chol_shock_cov = arma::trans(arma::chol(Sigma_hat));

        for (i=0; i<n_draws; i++) {
            beta_b  = beta_draws.slice(i);
            Sigma_b = Sigma_draws.slice(i);

            chol_shock_cov = arma::trans(arma::chol(Sigma_b));

            Y_forecast.zeros();
            Y_Th = Y_T;
            
            for (j=1; j<=horizon; j++) {
                Y_forecast.row(j-1) = Y_Th*beta_b + arma::trans(rmvnorm(chol_shock_cov,true));

                if (K_adj > M + c_int) {
                    Y_Th(0,arma::span(M+c_int,K_adj-1)) = Y_Th(0,arma::span(c_int,K_adj-M-1));
                }

                Y_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
            }
            //
            forecast_mat.slice(i) = Y_forecast;
        }
    } else {
        for (i=0; i<n_draws; i++) {
            beta_b = beta_draws.slice(i);

            Y_forecast.zeros();
            Y_Th = Y_T;

            for (j=1; j<=horizon; j++) {
                Y_forecast.row(j-1) = Y_Th*beta_b;

                if (K_adj > M + c_int) {
                    Y_Th(0,arma::span(M+c_int,K_adj-1)) = Y_Th(0,arma::span(c_int,K_adj-M-1));
                }

                Y_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
            }
            //
            forecast_mat.slice(i) = Y_forecast;
        }
    }
    //
    return forecast_mat;
}
