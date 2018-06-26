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
 * bvarm class
 */

#include "misc/misc.hpp"
namespace bm {
    #include "var/bvarm.hpp"
}

//
// data and basic setup

void
bm::bvarm::build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,nullptr,cons_term_inp,p_inp);
}

void
bm::bvarm::build(const arma::mat& data_raw, const arma::mat& data_ext, const bool cons_term_inp, const int p_inp)
{
    this->build_int(data_raw,&data_ext,cons_term_inp,p_inp);
}

void
bm::bvarm::build_int(const arma::mat& data_raw, const arma::mat* data_ext, const bool cons_term_inp, const int p_inp)
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
bm::bvarm::reset_draws()
{
    alpha_pt_mean.reset();
    alpha_pt_var.reset();

    beta_draws.reset();
}

//
// prior

void
bm::bvarm::prior(const arma::vec& coef_prior, double HP_1_inp, double HP_2_inp, double HP_3_inp)
{
    this->prior_int(coef_prior,nullptr,nullptr,HP_1_inp,HP_2_inp,HP_3_inp,nullptr);
}

void
bm::bvarm::prior(const arma::vec& coef_prior, double HP_1_inp, double HP_2_inp, double HP_3_inp, double HP_4_inp)
{
    this->prior_int(coef_prior,nullptr,nullptr,HP_1_inp,HP_2_inp,HP_3_inp,&HP_4_inp);
}

void
bm::bvarm::prior(const arma::vec& coef_prior, int var_type_inp, int decay_type_inp, 
                 double HP_1_inp, double HP_2_inp, double HP_3_inp, double HP_4_inp)
{
    this->prior_int(coef_prior,&var_type_inp,&decay_type_inp,HP_1_inp,HP_2_inp,HP_3_inp,&HP_4_inp);
}

void
bm::bvarm::prior_int(const arma::vec& coef_prior, const int* var_type_inp, const int* decay_type_inp, 
                     double HP_1_inp, double HP_2_inp, double HP_3_inp, const double* HP_4_inp)
{
    //
    // setup

    HP_1 = HP_1_inp;
    HP_2 = HP_2_inp;
    HP_3 = HP_3_inp;

    if (HP_4_inp) {
        HP_4 = *HP_4_inp;
    }

    var_type   = (var_type_inp) ? *var_type_inp : 1;
    decay_type = (decay_type_inp) ? *decay_type_inp : 1;

    //
    // OLS

    arma::mat beta_hat = arma::solve(X.t()*X,X.t()*Y);
    alpha_hat = arma::vectorise(beta_hat);

    arma::mat beta_pr_mean = arma::zeros(K,M);

    for (int i=0; i<M; i++) {
        beta_pr_mean(i + c_int,i) = coef_prior(i);
    }

    alpha_pr_mean = arma::vectorise(beta_pr_mean);

    //
    // residual variance

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
            X_AR.col(c_int + j) = X.col(c_int + i + j*M);
        }

        alpha_AR = arma::solve(X_AR.t()*X_AR,X_AR.t()*Y_AR);

        // ML estimate of the variance (no bias adjustment)

        arma::vec err_vec = Y_AR - X_AR*alpha_AR;
        sigma(i) = arma::dot(err_vec,err_vec) / static_cast<double>(n-p);
    }

    //
    // prior variance setup

    arma::mat beta_pr_var = arma::zeros(K,M);

    if (var_type==1) { // variance type 1
        if (cons_term) {
            beta_pr_var.row(0) = sigma.t() * HP_3;
        }

        for (int i=0; i < p; i++) {
            for (int j=0; j < M; j++) {
                for (int k=0; k < M; k++) {
                    beta_pr_var(i*M + k + c_int,j) = HP_2*sigma(j) / (std::pow(i+1,2)*sigma(k));
                }
                beta_pr_var(i*M + j + c_int,j) = HP_1 / std::pow(i+1,2);
            }
        }

        if (n_ext_vars > 0) {
            beta_pr_var.rows(K-n_ext_vars,K-1) = arma::trans(arma::repmat(sigma,1,n_ext_vars)) * HP_3;
        }
    } else if (var_type==2) { // variance type 2
        if (cons_term) {
            beta_pr_var.row(0).fill(HP_1 * HP_3);
        }

        double decay_i;

        for (int i=0; i < p; i++) {
            if (decay_type==1) {
                decay_i = decay_geo((double) i, HP_4);
            } else if (decay_type==2) {
                decay_i = decay_harm((double) i, HP_4);
            } else {
                printf("decay type not recognized.\n");
                return;
            }

            for (int j=0; j<M; j++) {
                for (int k=0; k<M; k++) {
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
        return;
    }

    //

    alpha_pr_var = arma::diagmat(arma::vectorise(beta_pr_var));
    Sigma_hat = arma::diagmat(sigma);
}

//
// posterior sampler

void
bm::bvarm::gibbs(const uint_t n_draws)
{
    beta_draws.set_size(K, M, n_draws);

    arma::mat inv_Sigma_hat = arma::inv(Sigma_hat);

    //

    arma::mat inv_alpha_pr_var = arma::inv(alpha_pr_var);

    alpha_pt_var = arma::inv_sympd(inv_alpha_pr_var + arma::kron(inv_Sigma_hat,X.t()*X));
    alpha_pt_mean = alpha_pt_var * (inv_alpha_pr_var*alpha_pr_mean + arma::vectorise(X.t() * Y * inv_Sigma_hat));

    arma::mat chol_alpha_pt_var = arma::chol(alpha_pt_var,"lower"); // lower triangular

    //

    if (only_stationary_draws)
    {
#ifdef BM_USE_OPENMP
        #pragma omp parallel for
#endif
        for (uint_t i=0; i < n_draws; i++)
        {
            arma::mat beta_draw = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean, chol_alpha_pt_var, true) ,K,M);

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
                    beta_draw = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean, chol_alpha_pt_var, true) ,K,M);
                }
            }

            beta_draws.slice(i) = beta_draw;
        }
    }
    else
    {
#ifdef BM_USE_OPENMP
        #pragma omp parallel for
#endif
        for (uint_t i=0; i < n_draws; i++) {
            beta_draws.slice(i) = arma::reshape( stats::rmvnorm<arma::mat>(alpha_pt_mean, chol_alpha_pt_var, true) ,K,M);
        }
    }
}

//
// IRFs

arma::cube
bm::bvarm::IRF(const uint_t n_irf_periods)
{
    return this->IRF_int(n_irf_periods, nullptr);
}

arma::cube
bm::bvarm::IRF(const uint_t n_irf_periods, const arma::mat& impact_mat)
{
    return this->IRF_int(n_irf_periods, &impact_mat);
}

arma::cube
bm::bvarm::IRF_int(const uint_t n_irf_periods, const arma::mat* impact_mat_inp)
{
    const uint_t n_draws = beta_draws.n_slices;
    const int K_adj = K - n_ext_vars;

    arma::mat impact_mat = (impact_mat_inp) ? *impact_mat_inp : arma::chol(Sigma_hat,"lower"); // impact_mat should be lower triangular

    arma::cube irfs(M, M, n_irf_periods*n_draws);

    //

    arma::mat Sigma_mat(M,M);
    if (irfs_lr_restrict) {
        Sigma_mat = (impact_mat_inp) ? impact_mat * impact_mat.t() : Sigma_hat;
    }

    arma::mat impact_mat_b(K_adj-c_int,M);
    arma::mat impact_mat_h(M,M);

#ifdef BM_USE_OPENMP
    #pragma omp parallel for firstprivate(impact_mat_b,impact_mat_h)
#endif
    for (uint_t j=1; j <= n_draws; j++)
    {
        arma::mat beta_b = beta_draws(arma::span(c_int,K_adj-1),arma::span(),arma::span(j-1,j-1));

        // long-run restrictions

        if (irfs_lr_restrict)
        {
            arma::mat poly_mat = arma::eye(M,M);
            for (int i=1; i<=p; i++) {
                poly_mat -= arma::trans(beta_b.rows(M*(i-1), M*i - 1));
            }

            arma::mat M_mat = arma::inv(poly_mat);

            impact_mat = poly_mat * arma::chol(M_mat*Sigma_mat*M_mat.t(),"lower");
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
bm::bvarm::FEVD(const uint_t n_periods)
{
    const uint_t n_draws = beta_draws.n_slices;
    const int K_adj = K - n_ext_vars;

    arma::mat impact_mat = arma::chol(Sigma_hat,"lower"); // impact_mat should be lower triangular

    arma::cube mse_cube(M, M, n_periods*n_draws);

    //

#ifdef BM_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint_t j=1; j <= n_draws; j++)
    {
        arma::mat beta_b = beta_draws(arma::span(c_int,K_adj-1),arma::span(),arma::span(j-1,j-1)); // b'th draw, minus coefficients on any external variables

        arma::mat poly_mat = arma::zeros(M,M);
        for (int i=1; i<=p; i++) {
            poly_mat += arma::trans(beta_b.rows(M*(i-1), M*i - 1));
        }

        if (irfs_lr_restrict)
        {   // long-run restrictions
            arma::mat M_mat = arma::inv(arma::eye(M,M) - poly_mat);

            impact_mat = (arma::eye(M,M) - poly_mat) * arma::chol(M_mat*Sigma_hat*M_mat.t(),"lower");
        }

        arma::mat mse_mat = Sigma_hat;
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
bm::bvarm::forecast(const uint_t horizon, const bool incl_shocks)
{
    return this->forecast_int(nullptr,horizon,incl_shocks);
}

arma::cube
bm::bvarm::forecast(const arma::mat& X_T, const uint_t horizon, const bool incl_shocks)
{
    return this->forecast_int(&X_T,horizon,incl_shocks);
}

arma::cube
bm::bvarm::forecast_int(const arma::mat* X_T_inp, const uint_t horizon, const bool incl_shocks)
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
    
    arma::cube forecast_cube(horizon, M, n_draws);

    //

    if (incl_shocks)
    {
        arma::mat chol_shock_cov = arma::chol(Sigma_hat,"lower");

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
                Y_forecast.row(j-1) = X_Th*beta_b + arma::trans(stats::rmvnorm<arma::mat>(arma::zeros(M,1),chol_shock_cov,true));

                if (K_adj > M + c_int) {
                    X_Th(0,arma::span(M+c_int,K_adj-1)) = X_Th(0,arma::span(c_int,K_adj-M-1));
                }

                X_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
            }

            //

            forecast_cube.slice(i) = Y_forecast;
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

            for (uint_t j=1; j<=horizon; j++)
            {
                Y_forecast.row(j-1) = X_Th*beta_b;

                if (K_adj > M + c_int) {
                    X_Th(0,arma::span(M+c_int,K_adj-1)) = X_Th(0,arma::span(c_int,K_adj-M-1));
                }

                X_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
            }

            //

            forecast_cube.slice(i) = Y_forecast;
        }
    }

    //

    return forecast_cube;
}
