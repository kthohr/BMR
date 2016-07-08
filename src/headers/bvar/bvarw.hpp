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
 * BVARW
 */

//#include "embed.hpp"
//#include "rmvnorm.hpp"
//#include "rinvwish.gpp"

// bvarw class
class bvarw
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
        arma::mat Z; // vec(Y) = Z alpha + vec(e)

        arma::mat alpha_pr_mean;  // prior mean
        arma::mat alpha_pr_var;   // prior variance

        arma::vec alpha_pt_mean;  // posterior mean
        arma::mat alpha_pt_var;   // posterior variance

        arma::mat Sigma_pr_scale; // prior scale matrix
        int Sigma_pr_dof;         // prior degrees of freedom

        arma::mat Sigma_pt_scale; // posterior scale matrix
        int Sigma_pt_dof;         // posterior degrees of freedom
        arma::mat Sigma_pt_mean;  // posterior mean

        arma::mat alpha_hat;      // OLS estimate of alpha
        arma::mat Sigma_hat;      // OLS-based estimation of covariance matrix of 'e'

        arma::cube beta_draws;    // posterior draws of beta
        arma::cube Sigma_draws;   // posterior draws of beta
        arma::cube irfs;          // irfs based on the posterior draws

        // member functions
        void data(arma::mat data_raw);
        void data(arma::mat data_raw, arma::mat data_ext);
        void prior(arma::vec coef_prior, double Xi_beta, arma::mat Xi_Sigma, int gamma);
        void prior(arma::vec coef_prior, arma::mat Xi_beta, arma::mat Xi_Sigma, int gamma);
        void gibbs(int n_draws, int n_burnin);
        void IRF(int n_irf_periods);
        arma::cube forecast(arma::mat Y_T, int horizon, bool incl_shocks);
        arma::cube forecast(arma::mat Y_T, int horizon, bool incl_shocks, arma::mat shock_cov);

    //private:
};

// This is broken into two functions (instead of using void data(arma::mat data_raw, arma::mat* data_ext))
// because Rcpp modules can't handle pointers
void bvarw::data(arma::mat data_raw)
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

void bvarw::data(arma::mat data_raw, arma::mat data_ext)
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

void bvarw::prior(arma::vec coef_prior, double Xi_beta, arma::mat Xi_Sigma, int gamma)
{
    int i;

    int np = n - p;
    //
    arma::mat beta_hat = arma::solve(X.t()*X,X.t()*Y);
    alpha_hat = arma::vectorise(beta_hat);

    arma::mat epsilon = Y - X*beta_hat;
    Sigma_hat = (epsilon.t() * epsilon) / ((double) np);
    //
    arma::mat beta_pr_mean = arma::zeros(K,M);

    for (i=0; i<M; i++) {
        beta_pr_mean(i + c_int,i) = coef_prior(i);
    }

    alpha_pr_mean = arma::vectorise(beta_pr_mean);
    //
    alpha_pr_var = arma::eye(K*M,K*M) * Xi_beta;
    //
    Sigma_pr_scale = Xi_Sigma;
    Sigma_pr_dof   = gamma;
    //
}

void bvarw::prior(arma::vec coef_prior, arma::mat Xi_beta, arma::mat Xi_Sigma, int gamma)
{
    int i;

    int np = n - p;
    //
    arma::mat beta_hat = arma::solve(X.t()*X,X.t()*Y);
    alpha_hat = arma::vectorise(beta_hat);

    arma::mat epsilon = Y - X*beta_hat;
    Sigma_hat = (epsilon.t() * epsilon) / ((double) np);
    //
    arma::mat beta_pr_mean = arma::zeros(K,M);

    for (i=0; i<M; i++) {
        beta_pr_mean(i + c_int,i) = coef_prior(i);
    }

    alpha_pr_mean = arma::vectorise(beta_pr_mean);
    //
    alpha_pr_var = Xi_beta;
    //
    Sigma_pr_scale = Xi_Sigma;
    Sigma_pr_dof   = gamma;
    //
}

void bvarw::gibbs(int n_draws, int n_burnin)
{
    int i;

    beta_draws.set_size(K, M, n_draws);
    Sigma_draws.set_size(M, M, n_draws);
    
    Sigma_pt_dof = n - p + Sigma_pr_dof;

    arma::mat XpX = X.t() * X;
    arma::mat XpY = X.t() * Y;
    arma::mat inv_alpha_pr_var = arma::inv(alpha_pr_var);
    //
    arma::mat inv_Sigma_b = arma::inv(Sigma_hat);
    arma::mat alpha_pt_var_b  = arma::inv_sympd(inv_alpha_pr_var + arma::kron(inv_Sigma_b,XpX));
    arma::vec alpha_pt_mean_b = alpha_pt_var_b * (inv_alpha_pr_var*alpha_pr_mean + arma::vectorise(XpY * inv_Sigma_b));
    
    arma::vec alpha_b = rmvnorm(alpha_pt_mean_b, alpha_pt_var_b, false);
    arma::mat beta_b = arma::reshape(alpha_b,K,M);
    //
    arma::mat epsilon = Y - X * beta_b;
    arma::mat Sigma_pt_scale_b = Sigma_pr_scale + epsilon.t() * epsilon;

    arma::mat Sigma_b = rinvwish(Sigma_pt_scale_b, Sigma_pt_dof);
    //
    for (i=0; i < (n_draws + n_burnin); i++) {
        inv_Sigma_b = arma::inv_sympd(Sigma_b);
        alpha_pt_var_b  = arma::inv_sympd(inv_alpha_pr_var + arma::kron(inv_Sigma_b,XpX));
        alpha_pt_mean_b = alpha_pt_var_b * (inv_alpha_pr_var*alpha_pr_mean + arma::vectorise(XpY * inv_Sigma_b));

        alpha_b = rmvnorm(alpha_pt_mean_b, alpha_pt_var_b, false);
        beta_b  = arma::reshape(alpha_b,K,M);
        //
        epsilon = Y - X * beta_b;
        Sigma_pt_scale_b = Sigma_pr_scale + epsilon.t() * epsilon;

        Sigma_b = rinvwish(Sigma_pt_scale_b, Sigma_pt_dof);
        //
        if (i >= n_burnin) {
            beta_draws.slice(i-n_burnin)  = beta_b;
            Sigma_draws.slice(i-n_burnin) = Sigma_b;
        }
    }
    //
    alpha_pt_mean = arma::vectorise(arma::mean(beta_draws,2));
    Sigma_pt_mean = arma::mean(Sigma_draws,2);
    //
}

void bvarw::IRF(int n_irf_periods)
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

arma::cube bvarw::forecast(arma::mat Y_T, int horizon, bool incl_shocks)
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

arma::cube bvarw::forecast(arma::mat Y_T, int horizon, bool incl_shocks, arma::mat shock_cov)
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
