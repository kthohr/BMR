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
 * BVARM
 */

//#include "embed.hpp"
//#include "rmvnorm.hpp"
#include "bvarm_aux.hpp"

// bvarm class
class bvarm
{
    public:
        bool cons_term; // if there is a constant (intercept) in the model

        int c_int;      // = 1 if cons_term == true
        int n;          // sample length (aka, 'T')
        int p;          // number of lags
        int M;          // number of endogenous variables
        int K;          // number of coefficients in each 
        int n_ext_vars; // number of 'external' variables

        int var_type;   // variance type (1 or 2)
        int decay_type; // decay type (1 or 2, for geometric and harmonic, resp.)

        arma::vec hyper_pars; // 'hyper' parameters

        arma::mat Y; // Y = X beta + e
        arma::mat X;
        arma::mat Z; // vec(Y) = Z alpha + vec(e)

        arma::mat alpha_pr_mean; // prior mean
        arma::mat alpha_pr_var;  // prior variance

        arma::vec alpha_pt_mean; // posterior mean
        arma::mat alpha_pt_var;  // posterior variance

        arma::mat alpha_hat;     // OLS estimate of alpha_hat
        arma::mat Sigma;         // covariance matrix of 'e' based on AR regressions

        arma::cube beta_draws;   // posterior draws of beta
        arma::cube irfs;         // irfs based on the beta draws

        // member functions
        void data(arma::mat data_raw);
        void data(arma::mat data_raw, arma::mat data_ext);
        int prior(arma::vec coef_prior);
        void gibbs(int n_draws);
        void IRF(int n_irf_periods);
        void IRF(arma::mat impact_mat, int n_irf_periods);
        arma::cube forecast(arma::mat Y_T, int horizon, bool incl_shocks);
        arma::cube forecast(arma::mat Y_T, int horizon, bool incl_shocks, arma::mat shock_cov);

    //private:
};

// This is broken into two functions (instead of using void data(arma::mat data_raw, arma::mat* data_ext))
// because Rcpp modules can't handle pointers
void bvarm::data(arma::mat data_raw)
{
    // data_raw is a n x M matrix with 'endogenous' variables
    // data_ext is a n x n_ext_vars matrix with 'exogenous' variables
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
    Z = arma::kron(arma::eye(M,M),X);
    //
}

void bvarm::data(arma::mat data_raw, arma::mat data_ext)
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
    Z = arma::kron(arma::eye(M,M),X);
    //
}

int bvarm::prior(arma::vec coef_prior)
{
    int i, j, k;

    int np = n-p;

    double HP_1 = hyper_pars(0);
    double HP_2 = hyper_pars(1);
    double HP_3 = hyper_pars(2);
    //
    arma::mat beta_hat = arma::solve(X.t()*X,X.t()*Y);
    alpha_hat = arma::vectorise(beta_hat);
    //
    arma::mat beta_pr_mean = arma::zeros(K,M);

    for (i=0; i<M; i++) {
        beta_pr_mean(i + c_int,i) = coef_prior(i);
    }

    alpha_pr_mean = arma::vectorise(beta_pr_mean);
    //
    // get 
    arma::vec sigma(M);

    arma::vec Y_AR(np), alpha_AR;
    arma::mat X_AR(np, c_int+p); 
    if (cons_term) {
        X_AR.col(0).fill(1);
    }

    for (i=0; i<M; i++) {
        Y_AR = Y.col(i);

        for (j=0; j<p; j++) {
            X_AR.col(c_int + j) = X.col(c_int + i + j*M);
        }
        //
        alpha_AR = arma::solve(X_AR.t()*X_AR,X_AR.t()*Y_AR);
        // ML estimate of the variance: (no bias adjustment)
        sigma(i) = arma::as_scalar(arma::trans(Y_AR - X_AR*alpha_AR) * (Y_AR - X_AR*alpha_AR)) / ((double) np);
    }
    //
    // prior variance setup
    arma::mat beta_pr_var = arma::zeros(K,M);

    if (var_type==1) {
        if (cons_term) {
            beta_pr_var.row(0) = sigma.t() * HP_3;
        }
        
        for (i=0; i<p; i++) {
            for (j=0; j<M; j++) {
                for (k=0; k<M; k++) {
                    beta_pr_var(i*M + k + c_int,j) = HP_2*sigma(j) / (std::pow(i+1,2)*sigma(k));
                }
                beta_pr_var(i*M + j + c_int,j) = HP_1 / std::pow(i+1,2);
            }
        }

        if (n_ext_vars > 0) {
            beta_pr_var.rows(K-n_ext_vars,K-1) = arma::trans(arma::repmat(sigma,1,n_ext_vars)) * HP_3;
        }
    } else if (var_type==2) {
        double HP_4 = hyper_pars(3);

        if (cons_term) {
            beta_pr_var.row(0).fill(HP_1 * HP_3);
        }

        double decay_i;
        
        for (i=0; i<p; i++) {
            if (decay_type==1) {
                decay_i = decay_geo((double) i, HP_4);
            } else if (decay_type==2) {
                decay_i = decay_harm((double) i, HP_4);
            } else {
                printf("decay type not recognized.\n");
                return -1;
            }

            for (j=0; j<M; j++) {
                for (k=0; k<M; k++) {
                    beta_pr_var(i*M + k + c_int,j) = HP_1*HP_2*sigma(k) / (decay_i*sigma(j));
                }
                beta_pr_var(i*M + j + c_int,j) = HP_1 / decay_i;
            }
        }

        if (n_ext_vars > 0) {
            beta_pr_var.rows(K-n_ext_vars,K-1).fill(HP_1 * HP_3);
        }
    } else {
        printf("variance type not recognized.\n");
        return -1;
    }
    //
    alpha_pr_var = arma::diagmat(arma::vectorise(beta_pr_var));
    Sigma = arma::diagmat(sigma);
    //
    return 0;
}

void bvarm::gibbs(int n_draws)
{
    int i;

    beta_draws.set_size(K, M, n_draws);

    arma::mat kron_inv_Sigma = arma::kron(arma::inv(Sigma),arma::eye(n-p,n-p));
    //
    arma::mat inv_alpha_pr_var = arma::inv(alpha_pr_var);
    
    alpha_pt_var = arma::inv_sympd(inv_alpha_pr_var + Z.t()*kron_inv_Sigma*Z);
    alpha_pt_mean = alpha_pt_var * (inv_alpha_pr_var*alpha_pr_mean + Z.t()*kron_inv_Sigma*arma::vectorise(Y));
    arma::mat chol_alpha_pt_var = arma::trans(arma::chol(alpha_pt_var)); // lower triangular

    arma::vec alpha(K*M);
    // no need for burnin
    for (i=0; i < n_draws; i++) {
        alpha = rmvnorm(alpha_pt_mean, chol_alpha_pt_var, true);

        beta_draws.slice(i) = arma::reshape(alpha,K,M);
    }
    //
}

void bvarm::IRF(int n_irf_periods)
{
    int i,j;

    int n_draws = beta_draws.n_slices;
    int K_adj = K - n_ext_vars;

    arma::mat impact_mat = arma::trans(arma::chol(Sigma));
    
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

void bvarm::IRF(arma::mat impact_mat, int n_irf_periods)
{
    /*
     * impact_mat should be lower triangular
     */
    int i,j;

    int n_draws = beta_draws.n_slices;
    int K_adj = K - n_ext_vars;
    
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

        impact_mat_b.zeros();
        impact_mat_b.rows(0,M-1) = impact_mat;

        irfs.slice((j-1)*n_irf_periods) = impact_mat;

        for (i=2; i<=n_irf_periods; i++) {
            impact_mat_h = beta_b_trans*impact_mat_b;
            irfs.slice((j-1)*n_irf_periods + (i-1)) = impact_mat_h;

            if (K_adj > M+c_int) {
                impact_mat_b.rows(M,K_adj-1-c_int) = impact_mat_b.rows(0,K_adj-M-1-c_int);
            }

            impact_mat_b.rows(0,M-1) = impact_mat_h;
        }
    }
    //
}

arma::cube bvarm::forecast(arma::mat Y_T, int horizon, bool incl_shocks)
{
    int i,j;

    int n_draws = beta_draws.n_slices;
    int K_adj = K - n_ext_vars;
    
    arma::mat beta_b(K_adj,M);       // bth draw
    arma::mat Y_forecast(horizon,M);
    arma::mat Y_Th = Y_T;

    arma::cube forecast_mat(horizon, M, n_draws);
    //
    if (incl_shocks) {
        arma::mat chol_shock_cov = arma::eye(M,M);

        for (i=0; i<n_draws; i++) {
            beta_b = beta_draws.slice(i);

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

arma::cube bvarm::forecast(arma::mat Y_T, int horizon, bool incl_shocks, arma::mat shock_cov)
{
    int i,j;

    int n_draws = beta_draws.n_slices;
    int K_adj = K - n_ext_vars;
    
    arma::mat beta_b(K_adj,M);       // bth draw
    arma::mat Y_forecast(horizon,M);
    arma::mat Y_Th = Y_T;

    arma::cube forecast_mat(horizon, M, n_draws);
    //
    if (incl_shocks) {
        arma::mat chol_shock_cov = arma::trans(arma::chol(shock_cov));

        for (i=0; i<n_draws; i++) {
            beta_b = beta_draws.slice(i);

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
