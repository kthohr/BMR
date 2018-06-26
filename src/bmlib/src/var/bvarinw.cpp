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
 * bvarinw class
 */

#include "misc/misc.hpp"
namespace bm {
    #include "var/bvarinw.hpp"
}

//
// data and basic setup

void
bm::bvarinw::build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,nullptr,cons_term_inp,p_inp);
}

void
bm::bvarinw::build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,&data_ext,cons_term_inp,p_inp);
}

void
bm::bvarinw::build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp)
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

    //

    arma::mat data_lagged = embed(data_raw,p);

    Y = data_lagged.cols(0,M-1);

    if (cons_term) {
        X = arma::join_rows(arma::ones(n-p),data_lagged.cols(M,data_lagged.n_cols-1));
    } else {
        X = data_lagged.cols(M,data_lagged.n_cols-1);
    }

    //

    if (data_ext && (n_ext_vars > 0)) {
        arma::mat X_ext = *data_ext;
        X_ext.shed_rows(0,p-1);

        X = arma::join_rows(X,X_ext);
    }

    // Z = arma::kron(arma::eye(M,M),X);
}

//
// reset

void
bm::bvarinw::reset_draws()
{
    alpha_pt_mean.reset();
    alpha_pt_var.reset();

    Sigma_pt_mean.reset();

    beta_draws.reset();
    Sigma_draws.reset();
}

//
// prior

void
bm::bvarinw::prior(const arma::vec& coef_prior, const double Xi_beta, const double Xi_Sigma, const int gamma)
{
    this->prior(coef_prior,Xi_beta*arma::eye(K*M,K*M),Xi_Sigma*arma::eye(M,M),gamma);
}

void
bm::bvarinw::prior(const arma::vec& coef_prior, const arma::mat& Xi_beta, const arma::mat& Xi_Sigma, const int gamma)
{
    //
    // MLE

    arma::mat beta_hat = arma::solve(X.t()*X,X.t()*Y);
    alpha_hat = arma::vectorise(beta_hat);

    const arma::mat epsilon = Y - X*beta_hat;
    Sigma_hat = (epsilon.t() * epsilon) / ((double) (n - p));

    //
    // prior mean of alpha

    arma::mat beta_pr_mean = arma::zeros(K,M);

    for (int i=0; i < M; i++) {
        beta_pr_mean(i + c_int,i) = coef_prior(i);
    }

    alpha_pr_mean = arma::vectorise(beta_pr_mean);

    alpha_pr_var = Xi_beta;

    //
    // error variance priors

    Sigma_pr_scale = Xi_Sigma;
    Sigma_pr_dof   = gamma;
}

//
// posterior sampler

void
bm::bvarinw::gibbs(const uint_t n_draws, const uint_t n_burnin)
{
    beta_draws.set_size(K, M, n_draws);
    Sigma_draws.set_size(M, M, n_draws);

    Sigma_pt_dof = n - p + Sigma_pr_dof;

    const arma::mat XpX = X.t() * X;
    const arma::mat XpY = X.t() * Y;

    const arma::mat inv_alpha_pr_var = arma::inv(alpha_pr_var);

    // first iteration

    arma::mat inv_Sigma_draw = arma::inv(Sigma_hat);

    arma::mat alpha_pt_var_b  = arma::inv_sympd(inv_alpha_pr_var + arma::kron(inv_Sigma_draw,XpX));
    arma::vec alpha_pt_mean_b = alpha_pt_var_b * (inv_alpha_pr_var*alpha_pr_mean + arma::vectorise(XpY * inv_Sigma_draw));

    arma::mat beta_draw = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean_b, alpha_pt_var_b), K,M);

    arma::mat epsilon = Y - X * beta_draw;
    arma::mat Sigma_draw = stats::rinvwish<arma::mat,double,int>(Sigma_pr_scale + epsilon.t()*epsilon, Sigma_pt_dof);

    //
    // begin loop

    for (uint_t i=0; i < (n_draws + n_burnin); i++)
    {
        inv_Sigma_draw = arma::inv_sympd(Sigma_draw);

        alpha_pt_var_b  = arma::inv_sympd(inv_alpha_pr_var + arma::kron(inv_Sigma_draw,XpX));
        alpha_pt_mean_b = alpha_pt_var_b * (inv_alpha_pr_var*alpha_pr_mean + arma::vectorise(XpY * inv_Sigma_draw));

        beta_draw = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean_b, alpha_pt_var_b), K,M);

        // check for a stationary draw

        if (only_stationary_draws)
        {
            bool loop_flag = true;

            while (loop_flag)
            {
                arma::mat poly_mat = arma::zeros(M,M);
                for (int i=1; i<=p; i++) {
                    poly_mat += arma::trans(beta_draw.rows(c_int + M*(i-1), c_int + M*i - 1));
                }

                arma::cx_vec eigvals = arma::eig_gen(poly_mat);

                if (arma::abs(eigvals).max() < 1.0) {
                    loop_flag = false; // escape
                }
                else {
                    beta_draw = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean_b, alpha_pt_var_b), K,M);
                }
            }
        }

        //

        epsilon = Y - X * beta_draw;
        Sigma_draw = stats::rinvwish<arma::mat,double,int>(Sigma_pr_scale + epsilon.t()*epsilon, Sigma_pt_dof);

        //

        if (i >= n_burnin)
        {
            beta_draws.slice(i-n_burnin)  = beta_draw;
            Sigma_draws.slice(i-n_burnin) = Sigma_draw;
        }
    }

    //

    alpha_pt_mean = arma::vectorise(arma::mean(beta_draws,2));
    Sigma_pt_mean = arma::mean(Sigma_draws,2);

    alpha_pt_var = arma::cov( cube_to_mat(beta_draws,true) );
}

//
// IRFs

arma::cube
bm::bvarinw::IRF(const uint_t n_irf_periods)
{
    const uint_t n_draws = beta_draws.n_slices;
    const int K_adj = K - n_ext_vars;

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

        // long-run restrictions

        if (irfs_lr_restrict)
        {
            arma::mat poly_mat = arma::eye(M,M);
            for (int i=1; i<=p; i++) {
                poly_mat -= arma::trans(beta_b.rows(M*(i-1), M*i - 1));
            }

            arma::mat M_mat = arma::inv(poly_mat);

            // arma::mat rotation_mat = arma::inv(impact_mat) * poly_mat * arma::chol(M_mat*Sigma_b*M_mat.t(),"lower");
            // impact_mat = impact_mat*rotation_mat;

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
bm::bvarinw::FEVD(const uint_t n_periods)
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
bm::bvarinw::forecast(const uint_t horizon, const bool incl_shocks)
{
    return this->forecast_int(nullptr,horizon,incl_shocks);
}

arma::cube
bm::bvarinw::forecast(const arma::mat& X_T, const uint_t horizon, const bool incl_shocks)
{
    return this->forecast_int(&X_T,horizon,incl_shocks);
}

arma::cube
bm::bvarinw::forecast_int(const arma::mat* X_T_inp, const uint_t horizon, const bool incl_shocks)
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
            arma::mat chol_shock_cov = arma::chol(Sigma_draws.slice(i),"lower");

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
