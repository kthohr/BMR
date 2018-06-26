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
 * bvars class
 */

#include "misc/misc.hpp"
namespace bm {
    #include "var/bvars.hpp"
}

//
// data and basic setup

void
bm::bvars::build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,nullptr,cons_term_inp,p_inp);
}

void
bm::bvars::build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,&data_ext,cons_term_inp,p_inp);
}

void
bm::bvars::build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp)
{
    // data_raw is a n x M matrix with 'endogenous' variables
    // data_ext is a n x n_ext_vars matrix with 'exogenous' variables

    cons_term = cons_term_inp;
    p = p_inp;

    c_int = !!cons_term;

    n = data_raw.n_rows;
    M = data_raw.n_cols;

    n_ext_vars = (data_ext) ? data_ext->n_cols : 0;

    K = M*p; // NOTE: K is defined differently to other BVAR models

    //

    arma::mat data_lagged = embed(data_raw,p);

    Y = data_lagged.cols(0,M-1);

    X = data_lagged.cols(M,data_lagged.n_cols-1);

    //

    arma::mat X_ext;

    if (data_ext && (n_ext_vars > 0)) {
        X_ext = *data_ext;
        q = n_ext_vars;

        if (cons_term) {
            q++;
            X_ext = arma::join_rows(arma::ones(n,1), X_ext);
        }
    } else {
        cons_term = true;
        q = 1;
        X_ext = arma::ones(n,1);
    }

    D = embed(X_ext,p);

    d   = D.cols(0,q-1);
    d_X = D.cols(q,q*(p+1)-1);

    D.cols(q,q*(p+1)-1) *= -1;
    //
}

//
// reset

void
bm::bvars::reset_draws()
{
    psi_pt_mean.reset();
    psi_pt_var.reset();

    alpha_pt_mean.reset();
    alpha_pt_var.reset();

    Sigma_pt_mean.reset();

    Psi_draws.reset();
    beta_draws.reset();
    Sigma_draws.reset();
}

//
// prior

arma::mat
bm::bvars::minn_pr_var()
{
    arma::vec sigma(M);

    arma::vec Y_AR(n-p), alpha_AR;
    arma::mat X_AR(n-p, c_int+p);

    if (cons_term) {
        X_AR.col(0).fill(1);
    }

    for (int i=0; i < M; i++)
    {
        Y_AR = Y.col(i);

        for (int j=0; j < p; j++) {
            X_AR.col(c_int + j) = X.col(i + j*M);
        }

        alpha_AR = arma::solve(X_AR.t()*X_AR,X_AR.t()*Y_AR);

        // ML estimate of the variance (no bias adjustment)

        arma::vec err_vec = Y_AR - X_AR*alpha_AR;
        sigma(i) = arma::dot(err_vec,err_vec) / static_cast<double>(n-p);
    }

    return arma::diagmat(sigma);
}

void
bm::bvars::prior(const arma::vec& coef_prior, const double HP_1, const double HP_2, 
                 const arma::mat& Psi_prior, const double Xi_psi, const int gamma,
                 const bool full_cov_prior)
{
    arma::mat Z = arma::join_rows(X,d);

    arma::mat var_coefs = arma::solve(Z.t()*Z,Z.t()*Y);

    arma::mat beta_hat = var_coefs.rows(0,K-1);
    arma::mat Phi_hat = var_coefs.rows(K,K+q-1);

    arma::mat poly_mat = arma::eye(M,M);
    for (int i=1; i<=p; i++) {
        poly_mat -= beta_hat.rows(M*(i-1), M*i - 1);
    }

    alpha_hat = arma::vectorise(beta_hat);
    psi_hat = arma::vectorise(Phi_hat*arma::inv(poly_mat));

    arma::mat epsilon = Y - Z*var_coefs;
    Sigma_hat = (epsilon.t() * epsilon) / ((double) (n - p));

    //

    psi_pr_mean = arma::vectorise(Psi_prior.t());
    psi_pr_var = Xi_psi*arma::eye(q*M,q*M);

    //

    arma::mat beta_pr_mean = arma::zeros(K,M);

    for (int i=0; i<M; i++) {
        beta_pr_mean(i,i) = coef_prior(i);
    }

    alpha_pr_mean = arma::vectorise(beta_pr_mean);

    //

    arma::mat Sigma_hat_pr = Sigma_hat;

    if (!full_cov_prior) {
        Sigma_hat_pr = minn_pr_var();
    }

    arma::mat beta_pr_var = arma::zeros(K,M);

    for (int i=0; i < p; i++) {
        for (int j=0; j < M; j++) {
            for (int k=0; k < M; k++) {
                beta_pr_var(i*M + k,j) = HP_2*Sigma_hat_pr(j,j) / (std::pow(i+1,2)*Sigma_hat_pr(k,k));
            }
            beta_pr_var(i*M + j,j) = HP_1 / std::pow(i+1,2);
        }
    }

    alpha_pr_var = arma::diagmat(arma::vectorise(beta_pr_var));

    //

    Sigma_pr_scale = Sigma_hat_pr;
    Sigma_pr_dof   = gamma;
}

//
// posterior sampler

void
bm::bvars::gibbs(const uint_t n_draws, const uint_t n_burnin)
{
    Psi_draws.set_size(q, M, n_draws);
    beta_draws.set_size(K, M, n_draws);
    Sigma_draws.set_size(M, M, n_draws);

    //

    Sigma_pt_dof = n - p + Sigma_pr_dof;

    arma::mat inv_psi_pr_var = arma::inv(psi_pr_var);

    arma::mat inv_alpha_pr_var = arma::inv(alpha_pr_var);
    arma::mat beta_b = arma::reshape(alpha_hat,K,M);

    arma::mat Sigma_b = Sigma_hat;
    arma::mat inv_Sigma_b = arma::inv_sympd(Sigma_b);

    //
    // Psi

    arma::mat beta_b_trans = beta_b.t();
    arma::mat U = arma::zeros(M*q + p*M*q,M*q);

    U.rows(0,M*q-1) = arma::eye(M*q,M*q);
    for (int j=1; j<=p; j++) {
        if (q > 1) {
            U.rows(M*q*j,M*q*(j+1)-1) = arma::kron(arma::eye(q,q),beta_b_trans.cols(M*(j-1),M*j-1));
        } else {
            U.rows(M*q*j,M*q*(j+1)-1) = beta_b_trans.cols(M*(j-1),M*j-1);
        }
    }

    arma::mat xi = Y - X*beta_b; // xi = PI(L) x_t

    arma::mat psi_pt_var_b  = arma::inv(U.t()*arma::kron(D.t()*D,inv_Sigma_b)*U + inv_psi_pr_var);
    arma::mat psi_pt_mean_b = psi_pt_var_b*(U.t()*arma::vectorise(inv_Sigma_b*xi.t()*D) + inv_psi_pr_var*psi_pr_mean);

    arma::mat Psi_b = arma::trans(arma::reshape( stats::rmvnorm<arma::mat>(psi_pt_mean_b, psi_pt_var_b), M,q));

    //
    // beta

    arma::mat Y_d = Y - d*Psi_b;
    arma::mat X_d = X - d_X*arma::kron(arma::eye(p,p),Psi_b);

    arma::mat alpha_pt_var_b  = arma::inv_sympd(inv_alpha_pr_var + arma::kron(inv_Sigma_b,X_d.t()*X_d));
    arma::vec alpha_pt_mean_b = alpha_pt_var_b * (inv_alpha_pr_var*alpha_pr_mean + arma::vectorise(X_d.t()*Y_d * inv_Sigma_b));

    beta_b = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean_b, alpha_pt_var_b), K,M);

    if (only_stationary_draws)
    {
        bool loop_flag = true;

        while (loop_flag)
        {
            arma::mat poly_mat = arma::zeros(M,M);
            for (int i=1; i<=p; i++) {
                poly_mat += arma::trans(beta_b.rows(M*(i-1), M*i - 1));
            }

            arma::cx_vec eigvals = arma::eig_gen(poly_mat);

            if (arma::abs(eigvals).max() < 1.0) {
                loop_flag = false; // escape
            }
            else {
                beta_b = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean_b, alpha_pt_var_b), K,M);
            }
        }
    }

    //
    // Sigma

    arma::mat epsilon = Y_d - X_d * beta_b;

    Sigma_b = stats::rinvwish<arma::mat,double,int>(Sigma_pr_scale + epsilon.t() * epsilon, Sigma_pt_dof);

    //
    // begin loop

    for (uint_t i=0; i < (n_draws + n_burnin); i++) 
    {
        beta_b_trans = beta_b.t();
        inv_Sigma_b = arma::inv_sympd(Sigma_b);

        U = arma::zeros(M*q + p*M*q,M*q);
        U.rows(0,M*q-1) = arma::eye(M*q,M*q);

        for (int j=1; j<=p; j++) {
            if (q > 1) {
                U.rows(M*q*j,M*q*(j+1)-1) = arma::kron(arma::eye(q,q),beta_b_trans.cols(M*(j-1),M*j-1));
            } else {
                U.rows(M*q*j,M*q*(j+1)-1) = beta_b_trans.cols(M*(j-1),M*j-1);
            }
        }

        xi = Y - X*beta_b; // PI(L) Y_t

        psi_pt_var_b  = arma::inv(U.t()*arma::kron(D.t()*D,inv_Sigma_b)*U + inv_psi_pr_var);
        psi_pt_mean_b = psi_pt_var_b*(U.t()*arma::vectorise(inv_Sigma_b*xi.t()*D) + inv_psi_pr_var*psi_pr_mean);

        Psi_b = arma::trans(arma::reshape( stats::rmvnorm<arma::mat>(psi_pt_mean_b, psi_pt_var_b), M,q));

        //

        Y_d = Y - d*Psi_b;
        X_d = X - d_X*arma::kron(arma::eye(p,p),Psi_b);

        alpha_pt_var_b  = arma::inv_sympd(inv_alpha_pr_var + arma::kron(inv_Sigma_b,X_d.t()*X_d));
        alpha_pt_mean_b = alpha_pt_var_b * (inv_alpha_pr_var*alpha_pr_mean + arma::vectorise(X_d.t()*Y_d * inv_Sigma_b));

        beta_b = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean_b, alpha_pt_var_b), K,M);

        //

        epsilon = Y_d - X_d * beta_b;

        Sigma_b = stats::rinvwish<arma::mat,double,int>(Sigma_pr_scale + epsilon.t() * epsilon, Sigma_pt_dof);

        //

        if (i >= n_burnin)
        {
            Psi_draws.slice(i-n_burnin)   = Psi_b;
            beta_draws.slice(i-n_burnin)  = beta_b;
            Sigma_draws.slice(i-n_burnin) = Sigma_b;
        }
    }

    //

    psi_pt_mean   = arma::vectorise(arma::mean(Psi_draws,2));
    alpha_pt_mean = arma::vectorise(arma::mean(beta_draws,2));
    Sigma_pt_mean = arma::mean(Sigma_draws,2);

    psi_pt_var   = arma::cov( cube_to_mat(Psi_draws,true) );
    alpha_pt_var = arma::cov( cube_to_mat(beta_draws,true) );

    //
}

//
// IRFs

arma::cube
bm::bvars::IRF(const uint_t n_irf_periods)
{
    const uint_t n_draws = beta_draws.n_slices;
    const int K_adj = K;

    arma::cube irfs(M, M, n_irf_periods*n_draws);

    //

    arma::mat impact_mat_b(K_adj,M);
    arma::mat impact_mat_h(M,M);

#ifdef BM_USE_OPENMP
    #pragma omp parallel for firstprivate(impact_mat_b,impact_mat_h)
#endif
    for (uint_t j=1; j <= n_draws; j++)
    {
        arma::mat beta_b = beta_draws(arma::span(0,K_adj-1),arma::span(),arma::span(j-1,j-1));

        arma::mat Sigma_b = Sigma_draws.slice(j-1);
        arma::mat impact_mat = arma::chol(Sigma_b,"lower");

        // long-run restrictions

        if (irfs_lr_restrict)
        {
            arma::mat poly_mat = arma::eye(M,M);
            for (int i=1; i<=p; i++) {
                poly_mat -= arma::trans(beta_b.rows(M*(i-1), M*i - 1));
            }

            arma::mat M_mat = arma::inv(poly_mat);
            impact_mat = poly_mat * arma::chol(M_mat*Sigma_b*M_mat.t(),"lower");
        }

        //

        impact_mat_b.zeros();
        impact_mat_b.rows(0,M-1) = impact_mat;

        irfs.slice((j-1)*n_irf_periods) = impact_mat;

        if (n_irf_periods > 1)
        {
            for (uint_t i=2; i <= n_irf_periods; i++)
            {
                impact_mat_h = beta_b.t()*impact_mat_b;
                irfs.slice((j-1)*n_irf_periods + (i-1)) = impact_mat_h;

                if (K_adj > M) {
                    impact_mat_b.rows(M,K_adj-1) = impact_mat_b.rows(0,K_adj-M-1);
                }

                impact_mat_b.rows(0,M-1) = std::move(impact_mat_h);
            }
        }
    }

    //

    return irfs;
}

//
// FEVD

arma::cube
bm::bvars::FEVD(const uint_t n_periods)
{
    const uint_t n_draws = beta_draws.n_slices;
    const int K_adj = K - n_ext_vars;

    arma::cube mse_cube(M, M, n_periods*n_draws);

    //

#ifdef BM_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint_t j=1; j <= n_draws; j++)
    {
        arma::mat beta_b = beta_draws(arma::span(0,K_adj-1),arma::span(),arma::span(j-1,j-1)); // b'th draw, minus coefficients on any external variables

        arma::mat Sigma_b = Sigma_draws.slice(j-1);
        arma::mat impact_mat = arma::chol(Sigma_b,"lower");

        arma::mat poly_mat = arma::zeros(M,M);
        for (int i=1; i<=p; i++) {
            poly_mat += arma::trans(beta_b.rows(M*(i-1), M*i - 1));
        }

        if (irfs_lr_restrict)
        {   // long-run restrictions
            arma::mat M_mat = arma::inv(arma::eye(M,M) - poly_mat);

            impact_mat = (arma::eye(M,M) - poly_mat) * arma::chol(M_mat*Sigma_b*M_mat.t(),"lower");
        }

        arma::mat mse_mat = Sigma_b;
        arma::mat fevd_mat = arma::pow(impact_mat,2);
        arma::mat iter_mat = impact_mat;

        arma::mat mse_slice = arma::zeros(M,M);

        for (int j=0; j < M; j++) {
            for (int k=0; k < M; k++) {
                mse_slice(j,k) = fevd_mat(j,k) / mse_mat(j,j);
            }
        }

        mse_cube.slice((j-1)*n_periods) = mse_slice;

        if (n_periods > 1)
        {
            for (uint_t i=2; i <= n_periods; i++)
            {
                iter_mat = poly_mat*iter_mat;

                mse_mat += iter_mat * iter_mat.t();

                for (int j=0; j < M; j++) {
                    for (int k=0; k < M; k++) {
                        fevd_mat(j,k) += std::pow(iter_mat(j,k),2);
                        mse_slice(j,k) = fevd_mat(j,k) / mse_mat(j,j);
                    }
                }

                mse_cube.slice((j-1)*n_periods + (i-1)) = mse_slice;
            }
        }
    }

    //

    return mse_cube;
}

//
// forecasting

arma::cube
bm::bvars::forecast(const uint_t horizon, const bool incl_shocks)
{
    return this->forecast_int(nullptr,horizon,incl_shocks);
}

arma::cube
bm::bvars::forecast(const arma::mat& X_T, const uint_t horizon, const bool incl_shocks)
{
    return this->forecast_int(&X_T,horizon,incl_shocks);
}

arma::cube
bm::bvars::forecast_int(const arma::mat* X_T_inp, const uint_t horizon, const bool incl_shocks)
{
    const uint_t n_draws = beta_draws.n_slices;
    const int K_adj = K;

    arma::mat X_T;

    if (X_T_inp) {
        X_T = *X_T_inp;
    } else {
        X_T = X.row(X.n_rows-1);

        if (K_adj > M) {
            X_T.cols(M,K_adj-1) = X_T.cols(0,K_adj-M-1);
        }

        X_T.cols(0,M-1) = Y.row(Y.n_rows-1);
    }

    arma::mat D_T = D.row(D.n_rows-1);

    arma::cube forecast_mat(horizon, M, n_draws);

    //

    if (incl_shocks)
    {
#ifdef BM_USE_OPENMP
        #pragma omp parallel for 
#endif
        for (uint_t i=0; i < n_draws; i++)
        {
            arma::mat Psi_b   = Psi_draws.slice(i);
            arma::mat beta_b  = beta_draws.slice(i);
            arma::mat Sigma_b = Sigma_draws.slice(i);

            // arma::mat D_term = D_T*arma::kron(arma::eye(q*(p+1),q*(p+1)),Psi_b)*arma::join_cols(arma::eye(M,M),beta_b); // Keith: check this
            arma::mat D_term = D_T*arma::kron(arma::eye(p+1,p+1),Psi_b)*arma::join_cols(arma::eye(M,M),beta_b); // Keith: check this

            arma::mat Y_forecast = arma::zeros(horizon,M);
            arma::mat X_Th = X_T;

            for (uint_t j=1; j <= horizon; j++)
            {
                Y_forecast.row(j-1) = X_Th*beta_b + D_term + arma::trans(stats::rmvnorm<arma::mat>(arma::zeros(M,1),Sigma_b));

                if (K_adj > M) {
                    X_Th(0,arma::span(M,K_adj-1)) = X_Th(0,arma::span(0,K_adj-M-1));
                }

                X_Th(0,arma::span(0,M-1)) = Y_forecast.row(j-1);
            }

            //

            forecast_mat.slice(i) = Y_forecast;
        }
    }
    else
    {
#ifdef BM_USE_OPENMP
        #pragma omp parallel for
#endif
        for (uint_t i=0; i < n_draws; i++)
        {
            arma::mat Psi_b  = Psi_draws.slice(i);
            arma::mat beta_b = beta_draws.slice(i);

            arma::mat D_term = D_T*arma::kron(arma::eye(p+1,p+1),Psi_b)*arma::join_cols(arma::eye(M,M),beta_b);

            arma::mat Y_forecast = arma::zeros(horizon,M);
            arma::mat X_Th = X_T;

            for (uint_t j=1; j<=horizon; j++)
            {
                Y_forecast.row(j-1) = X_Th*beta_b + D_term;

                if (K_adj > M) {
                    X_Th(0,arma::span(M,K_adj-1)) = X_Th(0,arma::span(0,K_adj-M-1));
                }

                X_Th(0,arma::span(0,M-1)) = Y_forecast.row(j-1);
            }

            //

            forecast_mat.slice(i) = Y_forecast;
        }
    }

    //

    return forecast_mat;
}
