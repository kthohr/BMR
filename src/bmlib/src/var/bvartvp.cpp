/*################################################################################
  ##
  ##   Copyright (C) 2011-2018 Keith O'Hara
  ##
  ##   This file is part of the BM++ C++ library.
  ##
  ##   BM++ is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   BM++ is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with BM++. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

/*
 * bvartvp class
 */

#include "misc/misc.hpp"
namespace bm {
    #include "var/bvartvp.hpp"
    #include "filter/dk_filter.hpp"
}

//
// data and basic setup

void
bm::bvartvp::build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,cons_term_inp,p_inp);
}

void
bm::bvartvp::build_int(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp)
{
    // data_raw is a n x M matrix with 'endogenous' variables
    // data_ext is a n x n_ext_vars matrix with 'exogenous' variables

    cons_term = cons_term_inp;
    p = p_inp;

    c_int = !!cons_term_inp;

    n = data_raw.n_rows;

    M = data_raw.n_cols;

    n_ext_vars = 0;

    K = c_int + M*p;
    //
    // arma::mat data_lagged = embed(data_raw,p);

    Y = data_raw; // not lagged yet; done in prior
}

//
// reset

void
bm::bvartvp::reset_draws()
{
    alpha_pt_mean.reset();
    Q_pt_mean.reset();
    Sigma_pt_mean.reset();

    alpha_draws.reset();
    Q_draws.reset();
    Sigma_draws.reset();
}

//
// prior

void
bm::bvartvp::prior(const int tau_inp, const double Xi_beta, const double Xi_Q_inp, const int gamma_Q, const double Xi_Sigma, const int gamma_S)
{
    tau = tau_inp;

    arma::mat Y_train = arma::trans(Y.rows(p,tau+p-1));

    // construct Z and Z_training

    Z = arma::zeros(M*(n-p),K*M);

    for (int i=(p+1); i <= n; i++) {
        int k = i - p - 1;

        if (cons_term) {
            Z(arma::span(k*M,(k+1)*M-1),arma::span(0,M-1)) = arma::eye(M,M);
        }

        for (int j=1; j <= p; j++) {
            Z(arma::span(k*M, (k+1)*M-1),arma::span((c_int-1)*M + M + M*M*(j-1), (c_int-1)*M + M + M*M*j - 1)) = arma::kron(arma::eye(M,M), Y.row(i-1-j));
        }
    }

    arma::mat Z_train = Z.rows(0,M*tau-1);
    Z.shed_rows(0,M*tau-1);

    //

    arma::mat alpha_denom = arma::zeros(K*M,K*M);
    arma::mat alpha_numer = arma::zeros(K*M,1);

    for (int i=0; i < tau; i++) {
        arma::mat Z_t = Z_train.rows(i*M,(i+1)*M-1);

        alpha_denom += Z_t.t()*Z_t;
        alpha_numer += Z_t.t()*Y_train.col(i);
    }

    alpha_pr_mean = arma::solve(alpha_denom,alpha_numer);
    alpha_hat = alpha_pr_mean;

    //

    arma::mat SSE = arma::zeros(M,M);

    for (int i=0; i < tau; i++) {
        arma::mat Z_t = Z_train.rows(i*M,(i+1)*M-1);

        SSE += (Y_train.col(i) - Z_t*alpha_pr_mean)*arma::trans(Y_train.col(i) - Z_t*alpha_pr_mean);
    }

    Sigma_hat = SSE / ((double) tau);

    //

    arma::mat inv_Sigma_hat = arma::inv(Sigma_hat);

    arma::mat inv_alpha_pr_var = arma::zeros(K*M,K*M);

    for (int i=0; i < tau; i++) {
        arma::mat Z_t = Z_train.rows(i*M,(i+1)*M-1);

        inv_alpha_pr_var += Z_t.t()*inv_Sigma_hat*Z_t;
    }

    alpha_pr_var = Xi_beta*arma::inv(inv_alpha_pr_var);

    //

    Xi_Q = Xi_Q_inp;

    Q_pr_scale = Xi_Q_inp*tau*arma::inv(inv_alpha_pr_var);
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

    arma::mat inv_alpha_pr_var = arma::inv(alpha_pr_var);

    // arma::mat Q_draw = Q_pr_scale / (double) Q_pr_dof;
    // arma::mat Q_chol = arma::chol(Q_draw,"lower");

    Q_pt_dof = n_adj + Q_pr_dof;
    Sigma_pt_dof = n_adj + Sigma_pr_dof;

    arma::mat Q_draw = Xi_Q*arma::eye(K*M,K*M);
    arma::mat Q_chol = std::sqrt(Xi_Q)*arma::eye(K*M,K*M);

    arma::mat Sigma_draw = Sigma_hat;
    arma::mat Sigma_chol = arma::chol(Sigma_draw,"lower");
    arma::mat inv_Sigma_draw = arma::inv(Sigma_draw);

    // storage

    alpha_draws  = arma::zeros(K*M,(n_adj+1),n_draws);
    Q_draws     = arma::zeros(K*M,K*M,n_draws);
    Sigma_draws = arma::zeros(M,M,n_draws);

    arma::mat alpha_tilde = arma::zeros(K*M,n_adj+1);

    //
    // alpha draw

    arma::mat alpha_pt_var = arma::zeros(K*M,K*M);
    arma::mat alpha_numer = arma::zeros(K*M,1);

    for (int i=0; i < n_adj; i++) {
        arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

        alpha_pt_var += Z_t.t() * inv_Sigma_draw * Z_t;
        alpha_numer += Z_t.t() * inv_Sigma_draw * (y.col(i) - Z_t*alpha_tilde.col(i));
    }

    alpha_pt_var = arma::inv_sympd(alpha_pt_var + inv_alpha_pr_var);
    arma::mat alpha_initial_draw = alpha_pt_var*(inv_alpha_pr_var*alpha_pr_mean + alpha_numer) + arma::chol(alpha_pt_var,"lower")*arma::randn(K*M,1);

    // run simulation smoother

    arma::mat epsilon_mat = arma::zeros(M,n_adj);
    arma::mat alpha_plus = arma::zeros(K*M,n_adj+1); // ensures that alpha_1 ~ N(0,Q)

    for (int i=0; i < n_adj; i++) {
        arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

        epsilon_mat.col(i) = y.col(i) - Z_t*alpha_initial_draw;
    }

    alpha_tilde = dk_filter(epsilon_mat, Z, Q_draw, Q_chol, Sigma_draw, Sigma_chol, n_adj, M, K);

    arma::mat alpha_draw = alpha_tilde + arma::repmat(alpha_initial_draw,1,n_adj+1);

    //
    // Q draw

    arma::mat Delta_alpha = arma::trans(alpha_draw.cols(1,n_adj)) - arma::trans(alpha_draw.cols(0,n_adj-1));

    arma::mat SSE_Q = arma::zeros(K*M,K*M);

    for (int i=0; i < n_adj; i++) {
        SSE_Q += Delta_alpha.row(i).t()*Delta_alpha.row(i);
    }

    arma::mat running_draw_mat = SSE_Q + Q_pr_scale;

    double scale_val = arma::abs(arma::diagvec(SSE_Q + Q_pr_scale)).min();

    running_draw_mat /= scale_val;

    Q_draw = stats::rinvwish<arma::mat>(running_draw_mat, Q_pt_dof);
    Q_draw = Q_draw * scale_val;
    Q_chol = arma::chol(Q_draw,"lower");

    //
    // Sigma draw

    epsilon_mat.zeros();
    arma::mat SSE_S = arma::zeros(M,M);

    for (int i=0; i < n_adj; i++) {
        epsilon_mat.col(i) = y.col(i) - Z.rows(i*M,(i+1)*M-1)*alpha_draw.col(i);

        SSE_S += epsilon_mat.col(i)*epsilon_mat.col(i).t();
    }

    running_draw_mat = SSE_S + Sigma_pr_scale;

    scale_val = arma::abs(arma::diagvec(running_draw_mat)).min();

    running_draw_mat /= scale_val;

    Sigma_draw = stats::rinvwish<arma::mat>(running_draw_mat, Sigma_pt_dof);

    Sigma_chol = std::sqrt(scale_val)*arma::chol(Sigma_draw,"lower");
    Sigma_draw = scale_val*Sigma_draw;

    //
    // begin loop
    //

    for (int rep=1; rep <= (n_draws + n_burnin); rep++) {

        inv_Sigma_draw = arma::inv(Sigma_draw);

        // alpha draw

        alpha_pt_var.zeros();
        alpha_numer.zeros();

        for (int i=0; i < n_adj; i++) {
            arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

            alpha_pt_var += Z_t.t() * inv_Sigma_draw * Z_t;
            alpha_numer += Z_t.t() * inv_Sigma_draw * (y.col(i) - Z_t*alpha_tilde.col(i));
        }

        alpha_pt_var = arma::inv_sympd(alpha_pt_var + inv_alpha_pr_var);

        alpha_initial_draw = alpha_pt_var*(inv_alpha_pr_var*alpha_pr_mean + alpha_numer) + arma::chol(alpha_pt_var,"lower")*arma::randn(K*M,1);

        // run simulation smoother

        epsilon_mat.zeros();

        alpha_plus.zeros();

        for (int i=0; i < n_adj; i++) {
            arma::mat Z_t = Z.rows(i*M,(i+1)*M-1);

            epsilon_mat.col(i) = y.col(i) - Z_t*alpha_initial_draw;
        }

        alpha_tilde = dk_filter(epsilon_mat, Z, Q_draw, Q_chol, Sigma_draw, Sigma_chol, n_adj, M, K);

        alpha_draw = alpha_tilde + arma::repmat(alpha_initial_draw,1,n_adj+1);

        //
        // Q draw

        Delta_alpha = arma::trans(alpha_draw.cols(1,n_adj)) - arma::trans(alpha_draw.cols(0,n_adj-1));

        SSE_Q = arma::zeros(K*M,K*M);
        for (int i=0; i < n_adj; i++) {
            SSE_Q += Delta_alpha.row(i).t()*Delta_alpha.row(i);
        }

        running_draw_mat = SSE_Q + Q_pr_scale;

        // scale_val = arma::abs(arma::diagvec(SSE_Q + Q_pr_scale)).min();
        scale_val = arma::nonzeros(arma::abs(running_draw_mat)).min(); // Q generally has conditioning issues, so we need to use a rescaling

        running_draw_mat /= scale_val;

        Q_draw = stats::rinvwish<arma::mat>(running_draw_mat, Q_pt_dof);
        Q_draw = Q_draw * scale_val;
        Q_draw = (Q_draw + Q_draw.t()) / 2.0;
        Q_chol = arma::chol(Q_draw,"lower");

        //
        // Sigma draw

        SSE_S = arma::zeros(M,M);

        for (int i=0; i < n_adj; i++) {
            epsilon_mat.col(i) = y.col(i) - Z.rows(i*M,(i+1)*M-1)*alpha_draw.col(i);

            SSE_S += epsilon_mat.col(i)*epsilon_mat.col(i).t();
        }

        running_draw_mat = SSE_S + Sigma_pr_scale;

        scale_val = arma::abs(arma::diagvec(running_draw_mat)).min();
        // scale_val = arma::abs(SSE_S + Sigma_pr_scale).min(); // bad idea

        running_draw_mat /= scale_val;

        Sigma_draw = stats::rinvwish<arma::mat>(running_draw_mat, Sigma_pt_dof);
        Sigma_draw = scale_val*Sigma_draw;
        Sigma_chol = arma::chol(Sigma_draw,"lower");

        // storage

        if (rep > n_burnin) {
            alpha_draws.slice(rep-n_burnin-1)  = alpha_draw;
            Q_draws.slice(rep-n_burnin-1)     = Q_draw;
            Sigma_draws.slice(rep-n_burnin-1) = Sigma_draw;
        }
    }
    //
    alpha_pt_mean = arma::mean(alpha_draws,2);
    Q_pt_mean = arma::mean(Q_draws,2);
    Sigma_pt_mean = arma::mean(Sigma_draws,2);
}

//
// IRFs

arma::cube
bm::bvartvp::IRF(const int n_irf_periods, const int time_ind)
{
    const int n_draws = alpha_draws.n_slices;
    const int K_adj = K - n_ext_vars;

    arma::cube irfs(M, M, n_irf_periods*n_draws);

    //

    arma::mat impact_mat_b(K_adj-c_int,M);
    arma::mat impact_mat_h(M,M);

    for (int j=1; j <= n_draws; j++) {
        arma::vec alpha_bt = alpha_draws(arma::span(),arma::span(time_ind,time_ind),arma::span(j-1,j-1));
        arma::mat beta_b = inside_trans(byrow(alpha_bt,c_int + p*M,M), cons_term); // put into standard BVAR 'beta' matrix format

        if (cons_term) {
            beta_b.shed_row(0);
        }

        arma::mat Sigma_b = Sigma_draws.slice(j-1);
        arma::mat impact_mat = arma::chol(Sigma_b,"lower");

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

            impact_mat_b.rows(0,M-1) = std::move(impact_mat_h);
        }
    }

    //

    return irfs;
}
