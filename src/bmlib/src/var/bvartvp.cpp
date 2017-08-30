/*################################################################################
  ##
  ##   Copyright (C) 2011-2017 Keith O'Hara
  ##
  ##   This file is part of the BMLib C++ library.
  ##
  ##   BMLib is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   BMLib is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ################################################################################*/

/*
 * bvartvp class
 */

#include "bmlib.hpp"

//
// data and basic setup

//
// prior

void
bm::bvartvp::build(const arma::mat& data_raw)
{
    this->build_int(data_raw);
}

// void
// bm::bvartvp::build(const arma::mat& data_raw, const arma::mat& data_ext)
// {
//     this->build_int(data_raw,&data_ext);
// }

void
bm::bvartvp::build_int(const arma::mat& data_raw)
{
    // data_raw is a n x M matrix with 'endogenous' variables
    // data_ext is a n x n_ext_vars matrix with 'exogenous' variables

    cons_term = true;

    c_int = 1;

    n = data_raw.n_rows;

    M = data_raw.n_cols;

    n_ext_vars = 0;

    K = c_int + M*p;
    //
    // arma::mat data_lagged = embed(data_raw,p);

    Y = data_raw; // not lagged yet; done in prior
}

void
bm::bvartvp::prior(const int tau_inp, const double Xi_beta, const double Xi_Q, const int gamma_Q, const double Xi_Sigma, const int gamma_S)
{
    tau = tau_inp;

    arma::mat Y_trsmp = arma::trans(Y.rows(p,tau+p-1));

    // construct Z and Z_training

    Z = arma::zeros(M*(n-p),K*M);

    for (int i=(p+1); i <= n; i++) {
        int k = i - p - 1;

        Z(arma::span(k*M,(k+1)*M-1),arma::span(0,M-1)) = arma::eye(M,M);

        for (int j=1; j<=p; j++) {
            Z(arma::span(k*M, (k+1)*M-1),arma::span(M+((M*M)*(j-1)), M-1+((M*M)*j))) = arma::kron(arma::eye(M,M), Y.row(i-1-j));
        }
    }

    arma::mat Z_trsmp = Z.rows(0,M*tau-1);
    Z.shed_rows(0,M*tau-1);

    //

    arma::mat beta_denom = arma::zeros(K*M,K*M);
    arma::mat beta_numer = arma::zeros(K*M,1);

    for (int i=0; i < tau; i++) {
        arma::mat Z_t = Z_trsmp.rows(i*M,(i+1)*M-1);

        beta_denom += Z_t.t()*Z_t;
        beta_numer += Z_t.t()*Y_trsmp.col(i);
    }

    beta_pr_mean = arma::solve(beta_denom,beta_numer);

    //

    arma::mat SSE = arma::zeros(M,M);

    for (int i=0; i < tau; i++) {
        arma::mat Z_t = Z_trsmp.rows(i*M,(i+1)*M-1);

        SSE += (Y_trsmp.col(i) - Z_t*beta_pr_mean)*arma::trans(Y_trsmp.col(i) - Z_t*beta_pr_mean);
    }

    Sigma_hat = SSE / ((double) tau);

    //

    arma::mat inv_Sigma_hat = arma::inv(Sigma_hat);

    arma::mat inv_beta_pr_var = arma::zeros(K*M,K*M);

    for (int i=0; i < tau; i++) {
        arma::mat Z_t = Z_trsmp.rows(i*M,(i+1)*M-1);

        inv_beta_pr_var += Z_t.t()*inv_Sigma_hat*Z_t;
    }

    beta_pr_var = Xi_beta*arma::inv(inv_beta_pr_var);

    //

    Q_pr_scale = Xi_Q*tau*arma::inv(inv_beta_pr_var);
    Q_pr_dof = gamma_Q;

    Sigma_pr_scale = arma::eye(M,M)*Xi_Sigma;
    Sigma_pr_dof = gamma_S;
}

void
bm::bvartvp::gibbs(const int n_draws, const int n_burnin)
{

    //
    // setup

    const int n_adj = n - tau - p;

    arma::mat y = arma::trans(Y.rows(tau+p,n-1));

    arma::mat inv_beta_pr_var = arma::inv(beta_pr_var);

    // arma::mat Q_draw = Q_pr_scale / (double) Q_pr_dof;
    // arma::mat Q_chol = arma::trans(arma::chol(Q_draw));

    // arma::mat Sigma_draw = Sigma_hat;
    // arma::mat Sigma_chol = arma::trans(arma::chol(Sigma_draw));
    // arma::mat inv_Sigma_draw = arma::inv(Sigma_draw);

    Q_pt_dof = n_adj + Q_pr_dof;
    Sigma_pt_dof = n_adj + Sigma_pr_dof;

    double consQ = 0.0001;
    arma::mat Q_draw = consQ*arma::eye(K*M,K*M);
    arma::mat Q_chol = std::sqrt(consQ)*arma::eye(K*M,K*M);
    arma::mat Sigma_draw = 0.1*arma::eye(M,M);
    arma::mat Sigma_chol = arma::trans(arma::chol(Sigma_draw));
    arma::mat inv_Sigma_draw = arma::inv(Sigma_draw);

    // storage

    beta_draws  = arma::zeros(K*M,(n_adj+1),n_draws);
    Q_draws     = arma::zeros(K*M,K*M,n_draws);
    Sigma_draws = arma::zeros(M,M,n_draws);

    arma::mat beta_tilde = arma::zeros(K*M,n_adj+1);

    //
    // beta draw

    std::cout << "begin gibbs" << std::endl;

    arma::mat beta_pt_var = arma::zeros(K*M,K*M);
    arma::mat beta_numer = arma::zeros(K*M,1);

    for (int i=0; i < n_adj; i++) {
        arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

        beta_pt_var += Z_t.t() * inv_Sigma_draw * Z_t;
        beta_numer += Z_t.t() * inv_Sigma_draw * (y.col(i) - Z_t*beta_tilde.col(i));
    }

    beta_pt_var = arma::inv_sympd(beta_pt_var + inv_beta_pr_var);
    arma::mat beta_initial_draw = beta_pt_var*(inv_beta_pr_var*beta_pr_mean + beta_numer) + arma::trans(arma::chol(beta_pt_var))*arma::randn(K*M,1);

    // run simulation smoother

    arma::mat epsilon_mat = arma::zeros(M,n_adj);
    arma::mat beta_plus = arma::zeros(K*M,n_adj+1); // ensures that beta_1 ~ N(0,Q)
    arma::mat Y_plus = arma::zeros(M,n_adj);

    for (int i=0; i < n_adj; i++) {
        arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

        epsilon_mat.col(i) = y.col(i) - Z_t*beta_initial_draw;
    }

    beta_tilde = dk_filter(epsilon_mat, Z, Q_draw, Q_chol, Sigma_draw, Sigma_chol, n_adj, M, K);

    arma::mat beta_draw = beta_tilde + arma::repmat(beta_initial_draw,1,n_adj+1);

    //
    // Q draw

    arma::mat Delta_beta = arma::trans(beta_draw.cols(1,n_adj)) - arma::trans(beta_draw.cols(0,n_adj-1));

    arma::mat SSE_Q = arma::zeros(K*M,K*M);

    for (int i=0; i < n_adj; i++) {
        SSE_Q += Delta_beta.row(i).t()*Delta_beta.row(i);
    }

    double scale_val = arma::abs(arma::diagvec(SSE_Q + Q_pr_scale)).min();

    Q_draw = stats::rinvwish((SSE_Q + Q_pr_scale) / scale_val, Q_pt_dof);
    Q_chol = std::sqrt(scale_val)*arma::trans(arma::chol(Q_draw));

    Q_draw = Q_draw * scale_val;

    //
    // Sigma draw

    epsilon_mat.zeros();
    arma::mat SSE_S = arma::zeros(M,M);

    for (int i=0; i < n_adj; i++) {
        epsilon_mat.col(i) = y.col(i) - Z.rows(i*M,(i+1)*M-1)*beta_draw.col(i);

        SSE_S += epsilon_mat.col(i)*epsilon_mat.col(i).t();
    }

    scale_val = arma::abs(arma::diagvec(SSE_S + Sigma_pr_scale)).min();

    Sigma_draw = stats::rinvwish((SSE_S + Sigma_pr_scale) / scale_val, Sigma_pt_dof);
    
    Sigma_chol = std::sqrt(scale_val)*arma::trans(arma::chol(Sigma_draw));
    Sigma_draw = scale_val*Sigma_draw;

    //
    // begin loop
    //

    for (int rep=1; rep <= (n_draws + n_burnin); rep++) {

        inv_Sigma_draw = arma::inv(Sigma_draw);

        // beta draw

        beta_pt_var.zeros();
        beta_numer.zeros();

        for (int i=0; i < n_adj; i++) {
            arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

            beta_pt_var += Z_t.t() * inv_Sigma_draw * Z_t;
            beta_numer += Z_t.t() * inv_Sigma_draw * (y.col(i) - Z_t*beta_tilde.col(i));
        }

        beta_pt_var = arma::inv_sympd(beta_pt_var + inv_beta_pr_var);

        beta_initial_draw = beta_pt_var*(inv_beta_pr_var*beta_pr_mean + beta_numer) + arma::trans(arma::chol(beta_pt_var))*arma::randn(K*M,1);

        // run simulation smoother

        epsilon_mat.zeros();

        beta_plus.zeros();
        Y_plus.zeros();

        for (int i=0; i < n_adj; i++) {
            arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

            epsilon_mat.col(i) = y.col(i) - Z_t*beta_initial_draw;
        }

        beta_tilde = dk_filter(epsilon_mat, Z, Q_draw, Q_chol, Sigma_draw, Sigma_chol, n_adj, M, K);

        beta_draw = beta_tilde + arma::repmat(beta_initial_draw,1,n_adj+1);

        //
        // Q draw

        Delta_beta = arma::trans(beta_draw.cols(1,n_adj)) - arma::trans(beta_draw.cols(0,n_adj-1));

        SSE_Q = arma::zeros(K*M,K*M);
        for (int i=0; i < n_adj; i++) {
            SSE_Q += Delta_beta.row(i).t()*Delta_beta.row(i);
        }

        // scale_val = arma::abs(arma::diagvec(SSE_Q + Q_pr_scale)).min();
        scale_val = arma::abs(SSE_Q + Q_pr_scale).min();

        Q_draw = stats::rinvwish((SSE_Q + Q_pr_scale) / scale_val, Q_pt_dof);
        Q_draw = Q_draw * scale_val;
        Q_draw = (Q_draw + Q_draw.t()) / 2.0;
        Q_chol = arma::trans(arma::chol(Q_draw));

        //
        // Sigma draw

        SSE_S = arma::zeros(M,M);

        for (int i=0; i < n_adj; i++) {
            epsilon_mat.col(i) = y.col(i) - Z.rows(i*M,(i+1)*M-1)*beta_draw.col(i);

            SSE_S += epsilon_mat.col(i)*epsilon_mat.col(i).t();
        }

        // scale_val = arma::abs(arma::diagvec(SSE_S + Sigma_pr_scale)).min();
        scale_val = arma::abs(SSE_S + Sigma_pr_scale).min();

        Sigma_draw = stats::rinvwish((SSE_S + Sigma_pr_scale) / scale_val, Sigma_pt_dof);
        Sigma_draw = scale_val*Sigma_draw;
        Sigma_chol = arma::trans(arma::chol(Sigma_draw));

        // storage

        if (rep > n_burnin) {
            beta_draws.slice(rep-n_burnin-1)  = beta_draw;
            Q_draws.slice(rep-n_burnin-1)     = Q_draw;
            Sigma_draws.slice(rep-n_burnin-1) = Sigma_draw;
        }
    }
    //
    beta_pt_mean = arma::mean(beta_draws,2);
    Q_pt_mean = arma::mean(Q_draws,2);
    Sigma_pt_mean = arma::mean(Sigma_draws,2);
}

//
// IRFs

void
bm::bvartvp::IRF(const int n_irf_periods)
{
    const int n_draws = beta_draws.n_slices;
    const int K_adj = K - n_ext_vars;
    
    irfs.set_size(M, M, n_irf_periods*n_draws);
    //
    arma::mat impact_mat_b(K_adj-c_int,M);
    arma::mat impact_mat_h(M,M);
    
    for (int j=1; j <= n_draws; j++) {
        arma::mat beta_b = beta_draws(arma::span(c_int,K_adj-1),arma::span(),arma::span(j-1,j-1)); // b'th draw, minus coefficients on any external variables 

        arma::mat Sigma_b = Sigma_draws.slice(j-1);
        arma::mat impact_mat = arma::trans(arma::chol(Sigma_b));
        //
        impact_mat_b.zeros();
        impact_mat_b.rows(0,M-1) = impact_mat;

        irfs.slice((j-1)*n_irf_periods) = impact_mat;

        for (int i=2; i <= n_irf_periods; i++) {
            impact_mat_h = beta_b.t()*impact_mat_b;
            irfs.slice((j-1)*n_irf_periods + (i-1)) = impact_mat_h;

            if(K_adj > M + c_int){
                impact_mat_b.rows(M,K_adj-c_int-1) = impact_mat_b.rows(0,K_adj-M-c_int-1);
            }

            impact_mat_b.rows(0,M-1) = impact_mat_h;
        }
    }
    //
}
