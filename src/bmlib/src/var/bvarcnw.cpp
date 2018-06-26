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
 * bvarcnw class
 */

#include "misc/misc.hpp"
namespace bm {
    #include "var/bvarcnw.hpp"
}

//
// data and basic setup

void
bm::bvarcnw::build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,nullptr,cons_term_inp,p_inp);
}

void
bm::bvarcnw::build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,&data_ext,cons_term_inp,p_inp);
}

void
bm::bvarcnw::build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp)
{
    // data_raw is a n x M matrix with 'endogenous' variables
    // data_ext is a n x n_ext_vars matrix with 'exogenous' variables

    cons_term = cons_term_inp;
    p = p_inp;

    c_int = !!cons_term;

    n = data_raw.n_rows;
    M = data_raw.n_cols;

    n_ext_vars = (data_ext) ? data_ext->n_cols : 0;

    K = c_int + M*p + n_ext_vars;
    K_adj = c_int + M*p;

    //

    arma::mat data_lagged = embed(data_raw,p);

    Y = data_lagged.cols(0,M-1);

    if (cons_term) {
        X = arma::join_rows(arma::ones(n-p),data_lagged.cols(M,data_lagged.n_cols-1));
    } else {
        X = data_lagged.cols(M,data_lagged.n_cols-1);
    }

    //

    if (data_ext && (n_ext_vars > 0))
    {
        arma::mat X_ext = *data_ext;
        X_ext.shed_rows(0,p-1);

        X = arma::join_rows(X,X_ext);
    }
}

//
// reset

void
bm::bvarcnw::reset_draws()
{
    alpha_pt_mean.reset();
    alpha_pt_var.reset();

    Sigma_pt_mean.reset();

    beta_draws.reset();
    Sigma_draws.reset();
}

//
// prior

arma::mat
bm::bvarcnw::minn_pr_var()
{
    if (p_AR < 0) {
        p_AR = p;
    }

    arma::vec sigma(M);

    arma::vec Y_AR(n-p), alpha_AR;
    arma::mat X_AR(n-p, c_int+p_AR);

    if (cons_term) {
        X_AR.col(0).fill(1);
    }

    for (int i=0; i < M; i++)
    {
        arma::vec Y_AR = Y.col(i);

        for (int j=0; j < p_AR; j++) {
            X_AR.col(c_int + j) = X.col(i + j*M);
        }

        alpha_AR = arma::solve(X_AR.t()*X_AR,X_AR.t()*Y_AR);

        // ML estimate of the variance (no bias adjustment)

        arma::vec err_vec = Y_AR - X_AR*alpha_AR;
        sigma(i) = arma::dot(err_vec,err_vec) / static_cast<double>(n-p_AR);
    }

    return arma::diagmat(sigma);
}

void
bm::bvarcnw::prior(const arma::vec& coef_prior, const int gamma)
{
    this->prior(coef_prior,HP_1,HP_3,gamma);
}

void
bm::bvarcnw::prior(const arma::vec& coef_prior, double HP_1_inp, double HP_3_inp, const int gamma,
                   const bool full_cov_prior)
{
    //
    // MLE

    arma::mat beta_hat = arma::solve(X.t()*X,X.t()*Y);
    alpha_hat = arma::vectorise(beta_hat);

    const arma::mat epsilon = Y - X*beta_hat;
    Sigma_hat = (epsilon.t() * epsilon) / static_cast<double> (n - p);

    //
    // prior mean of alpha

    beta_pr_mean = arma::zeros(K,M);

    for (int i=0; i < M; i++) {
        beta_pr_mean(i + c_int,i) = coef_prior(i);
    }

    alpha_pr_mean = arma::vectorise(beta_pr_mean);

    //
    //

    arma::mat Sigma_hat_pr = Sigma_hat;

    if (!full_cov_prior) {
        Sigma_hat_pr = minn_pr_var();
    }

    HP_1 = HP_1_inp;
    HP_3 = HP_3_inp;

    beta_pr_var = arma::zeros(K,K); // variance assumed the same for each equation

    if (cons_term) {
        beta_pr_var(0,0) = HP_1 * HP_3;
    }

    for (int i=0; i < p; i++) {
        for (int j=0; j < M; j++) {
            beta_pr_var(c_int + i*M + j, c_int + i*M + j) = HP_1 / ( Sigma_hat_pr(j,j) * std::pow(i+1,2) );
        }
    }

    if (n_ext_vars > 0) {
        for (int j=K_adj; j < K; j++) {
            beta_pr_var(j,j) = HP_1 * HP_3;
        }
    }

    alpha_pr_var = beta_pr_var;

    //
    // error variance priors

    Sigma_pr_scale = Sigma_hat_pr;
    Sigma_pr_dof   = gamma;
}

//
// posterior sampler

void
bm::bvarcnw::gibbs(const uint_t n_draws)
{
    beta_draws.set_size(K, M, n_draws);
    Sigma_draws.set_size(M, M, n_draws);

    Sigma_pt_dof = n - p + Sigma_pr_dof;

    const arma::mat XpX = X.t() * X;
    const arma::mat XpY = X.t() * Y;

    // first iteration

    arma::mat inv_beta_pr_var = arma::inv(beta_pr_var);

    arma::mat Qa_hat = inv_beta_pr_var + XpX;
    arma::mat invQa_hat = arma::inv(Qa_hat);
    arma::mat chol_invQa_hat = arma::chol(invQa_hat,"lower");

    arma::mat beta_pt_mean = invQa_hat * (inv_beta_pr_var*beta_pr_mean + XpY);
    arma::mat alpha_hat = arma::vectorise(beta_pt_mean);

    arma::mat Qe_hat_adj = arma::inv( Sigma_pr_dof*Sigma_pr_scale + beta_pr_mean.t()*inv_beta_pr_var*beta_pr_mean \
                                        + Y.t()*Y - beta_pt_mean.t() * Qa_hat * beta_pt_mean );

    arma::mat chol_Qe_hat_adj = arma::chol(Qe_hat_adj,"lower");

    //
    // begin loop

#ifdef BM_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint_t i=0; i < n_draws; i++)
    {
        Sigma_draws.slice(i) = arma::inv(stats::rwish<arma::mat>(chol_Qe_hat_adj,Sigma_pt_dof,true));
    }

#ifdef BM_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint_t i=0; i < n_draws; i++)
    {
        arma::mat chol_Sigma = arma::chol(Sigma_draws.slice(i),"lower");
        arma::mat chol_alpha_pt_var = arma::kron(chol_Sigma,chol_invQa_hat); // kronecker of chol = chol of kronecker

        arma::mat beta_draw = arma::reshape( stats::rmvnorm<arma::mat>(alpha_hat, chol_alpha_pt_var, true), K,M);

        //

        if (only_stationary_draws)
        {
            bool loop_flag = true;

            while (loop_flag)
            {
                // arma::mat poly_mat = arma::zeros(M,M);
                // for (int i=1; i<=p; i++) {
                //     poly_mat += arma::trans(beta_draw.rows(c_int + M*(i-1), c_int + M*i - 1));
                // }

                // arma::cx_vec eigvals = arma::eig_gen(poly_mat);

                arma::cx_vec eigvals = arma::eig_gen(companion_form_matrix(beta_draw,c_int,K_adj));

                if (arma::abs(eigvals).max() < 1.0) {
                    loop_flag = false; // escape
                }
                else {
                    beta_draw = arma::reshape( stats::rmvnorm<arma::mat>(alpha_hat, chol_alpha_pt_var, true), K,M);
                }
            }
        }

        //

        beta_draws.slice(i) = beta_draw;
    }

    //

    alpha_pt_mean = arma::vectorise(arma::mean(beta_draws,2));
    Sigma_pt_mean = arma::mean(Sigma_draws,2);

    alpha_pt_var = arma::cov( cube_to_mat(beta_draws,true) );
}

//
// IRFs

arma::cube
bm::bvarcnw::IRF(const uint_t n_irf_periods)
{
    const uint_t n_draws = beta_draws.n_slices;

    arma::cube irfs(M, M, n_irf_periods*n_draws);

    //

    arma::mat impact_mat_b(K_adj-c_int,M);
    arma::mat impact_mat_h(M,M);

#ifdef BM_USE_OPENMP
    #pragma omp parallel for firstprivate(impact_mat_b,impact_mat_h)
#endif
    for (uint_t j=1; j <= n_draws; j++)
    {
        arma::mat beta_b = beta_draws(arma::span(c_int,K_adj-1),arma::span(),arma::span(j-1,j-1)); // b'th draw, minus coefficients on any external variables

        arma::mat Sigma_b = Sigma_draws.slice(j-1);
        arma::mat impact_mat = arma::chol(Sigma_b,"lower");

        if (irfs_lr_restrict)
        {   // long-run restrictions
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

                if (K_adj > M + c_int) {
                    impact_mat_b.rows(M,K_adj-c_int-1) = impact_mat_b.rows(0,K_adj-M-c_int-1);
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
bm::bvarcnw::FEVD(const uint_t n_periods)
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
        arma::mat beta_b = beta_draws(arma::span(c_int,K_adj-1),arma::span(),arma::span(j-1,j-1)); // b'th draw, minus coefficients on any external variables

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
bm::bvarcnw::forecast(const uint_t horizon, const bool incl_shocks)
{
    return this->forecast_int(nullptr,horizon,incl_shocks);
}

arma::cube
bm::bvarcnw::forecast(const arma::mat& X_T, const uint_t horizon, const bool incl_shocks)
{
    return this->forecast_int(&X_T,horizon,incl_shocks);
}

arma::cube
bm::bvarcnw::forecast_int(const arma::mat* X_T_inp, const uint_t horizon, const bool incl_shocks)
{
    const uint_t n_draws = beta_draws.n_slices;
    const int K_adj = K - n_ext_vars;

    arma::mat X_T;

    if (X_T_inp) {
        X_T = *X_T_inp;
    } else {
        X_T = X.row(X.n_rows-1);

        if (K_adj > M + c_int) {
            X_T.cols(c_int+M,K_adj-1) = X_T.cols(c_int,K_adj-M-1);
        }

        X_T.cols(c_int,c_int+M-1) = Y.row(Y.n_rows-1);
    }

    arma::cube forecast_mat(horizon, M, n_draws);
    
    //

    if (incl_shocks)
    {
#ifdef BM_USE_OPENMP
        #pragma omp parallel for 
#endif
        for (uint_t i=0; i < n_draws; i++)
        {
            arma::mat beta_b  = beta_draws.slice(i);
            arma::mat Sigma_b = Sigma_draws.slice(i);

            arma::mat chol_shock_cov = arma::chol(Sigma_b,"lower");

            arma::mat Y_forecast = arma::zeros(horizon,M);
            arma::mat X_Th = X_T;

            for (uint_t j=1; j <= horizon; j++)
            {
                Y_forecast.row(j-1) = X_Th*beta_b + arma::trans(stats::rmvnorm<arma::mat>(arma::zeros(M,1),chol_shock_cov,true));

                if (K_adj > M + c_int) {
                    X_Th(0,arma::span(M+c_int,K_adj-1)) = X_Th(0,arma::span(c_int,K_adj-M-1));
                }

                X_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
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
            arma::mat beta_b = beta_draws.slice(i);

            arma::mat Y_forecast = arma::zeros(horizon,M);
            arma::mat X_Th = X_T;

            for (uint_t j=1; j <= horizon; j++)
            {
                Y_forecast.row(j-1) = X_Th*beta_b;

                if (K_adj > M + c_int) {
                    X_Th(0,arma::span(M+c_int,K_adj-1)) = X_Th(0,arma::span(c_int,K_adj-M-1));
                }

                X_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
            }

            //
            
            forecast_mat.slice(i) = Y_forecast;
        }
    }

    //

    return forecast_mat;
}
