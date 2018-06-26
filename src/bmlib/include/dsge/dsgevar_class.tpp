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
 * DSGE-VAR model class
 */

template<typename T>
void
dsgevar<T>::build(const arma::mat& data_raw, const bool cons_term_inp, const int p_inp, const double lambda_inp)
{
    // data_raw is a n x M matrix

    cons_term = cons_term_inp;

    p = p_inp;

    c_int = !!cons_term;

    n = data_raw.n_rows;
    M = data_raw.n_cols;

    K = c_int + M*p;

    lambda = (std::isfinite(lambda_inp)) ? std::max(lambda_inp, (double) (M*(p+1) + c_int) / (double) n) : lambda_inp;

    //

    dsge_obj.estim_data = data_raw;

    arma::mat data_lagged = embed(data_raw,p);

    Y = data_lagged.cols(0,M-1);

    if (cons_term) {
        X = arma::join_rows(arma::ones(n-p),data_lagged.cols(M,data_lagged.n_cols-1));
    } else {
        X = data_lagged.cols(M,data_lagged.n_cols-1);
    }

    //

    YY = Y.t()*Y / (double) (n-p);
    XY = X.t()*Y / (double) (n-p);
    XX = X.t()*X / (double) (n-p);
}

//
// set bounds and priors

template<typename T>
void
dsgevar<T>::set_bounds(const arma::vec& lower_bounds_inp, const arma::vec& upper_bounds_inp)
{
    dsge_obj.lower_bounds = lower_bounds_inp;
    dsge_obj.upper_bounds = upper_bounds_inp;
}

template<typename T>
void
dsgevar<T>::set_prior(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp)
{
    dsge_obj.set_prior(prior_form_inp,prior_pars_inp);
}

//
// model-implied moments

template<typename T>
void
dsgevar<T>::model_moments(arma::mat& Gamma_YY, arma::mat& Gamma_YX, arma::mat& Gamma_XX)
const
{
    this->model_moments(dsge_obj, Gamma_YY, Gamma_YX, Gamma_XX);
}

template<typename T>
void
dsgevar<T>::model_moments(const dsge<T>& dsge_model_inp, arma::mat& Gamma_YY, arma::mat& Gamma_YX, arma::mat& Gamma_XX)
const
{
    // we explicitly include the DSGE object to allow for parallelized gibbs sampling
    arma::mat KMC2 = dsge_model_inp.kalman_mat_C*dsge_model_inp.kalman_mat_C.t();

    const arma::mat Sigma_SS = lyapunov_dbl(dsge_model_inp.kalman_mat_Q, dsge_model_inp.kalman_mat_F);

    //

    arma::cube Sigma_X = arma::zeros(M, M, p+1);

    Sigma_X.slice(0) = KMC2 + dsge_model_inp.kalman_mat_H.t()*Sigma_SS*dsge_model_inp.kalman_mat_H + dsge_model_inp.kalman_mat_R;

    arma::mat Fp = arma::eye(dsge_model_inp.kalman_mat_F.n_cols,dsge_model_inp.kalman_mat_F.n_cols);

    for (int i=1; i<=p; i++) {
        Fp *= dsge_model_inp.kalman_mat_F;
        Sigma_X.slice(i) = KMC2 + dsge_model_inp.kalman_mat_H.t()*Fp*Sigma_SS*dsge_model_inp.kalman_mat_H;
    }

    //
    // model-implied moments

    Gamma_YY = Sigma_X.slice(0);

    Gamma_XX = arma::zeros(M*p,M*p);
    Gamma_YX = arma::zeros(M,M*p);
    arma::mat Gamma_CY = arma::zeros(1,M*p);

    arma::uvec r_inds, c_inds;

    for (int i=0; i < p; i++) {
        r_inds = uvec_linspace(i*M,(i+1)*M - 1);

        for (int j=0; j < p; j++) {
            c_inds = uvec_linspace(j*M,(j+1)*M - 1);

            if (i == j) {
                Gamma_XX(r_inds,c_inds) = Sigma_X.slice(0);
            } else if (i < j) {
                Gamma_XX(r_inds,c_inds) = Sigma_X.slice(j-i);
            } else {
                Gamma_XX(r_inds,c_inds) = Sigma_X.slice(i-j);
            }
        }

        Gamma_CY.cols(r_inds) = dsge_model_inp.kalman_mat_C.t();
        Gamma_YX.cols(r_inds) = Sigma_X.slice(i+1);
    }

    if (cons_term) {
        Gamma_YX = arma::join_rows(dsge_model_inp.kalman_mat_C,Gamma_YX);
        Gamma_XX = arma::join_cols( arma::join_rows(arma::ones(Gamma_CY.n_rows,1),Gamma_CY), arma::join_rows(Gamma_CY.t(), Gamma_XX) );
    }
}

//
// log likelihood

template<typename T>
double
dsgevar<T>::log_likelihood(const arma::mat& Gamma_YY, const arma::mat& Gamma_YX, const arma::mat& Gamma_XX, const arma::mat& Gamma_bar_YY, const arma::mat& Gamma_bar_YX, const arma::mat& Gamma_bar_XX)
const
{
    // note: assumes finite lambda value
    const double lambdaT = (n-p)*lambda;
    const double lambda_recip = 1.0/lambda;

    const arma::mat Sigma_eps = Gamma_YY - Gamma_YX*arma::inv_sympd(Gamma_XX)*Gamma_YX.t();
    const arma::mat Sigma_hat_eps = Gamma_bar_YY - Gamma_bar_YX.t()*arma::inv_sympd(Gamma_bar_XX)*Gamma_bar_YX;

    //

    const double cons_val = - M*(n-p)*std::log(lambdaT*arma::datum::pi)/2.0 + log_GPR(M, (1.0+lambda)*(n-p) - M*p - c_int, lambda*(n-p) - M*p - c_int);

    const double term_1 = M * ( - std::log(arma::det(Gamma_XX + lambda_recip*XX)) + std::log(arma::det(Gamma_XX)) ) / 2.0;
    // const double term_2 = - (n - p + lambdaT - M*p - c_int)*std::log(arma::det((1.0 + lambda_recip)*Sigma_hat_eps)) / 2.0;
    const double term_2 = - (n - p + lambdaT - M*p - c_int)*( M*std::log(1.0 + lambda_recip) + std::log(arma::det(Sigma_hat_eps)) ) / 2.0; // avoids bug in BMR
    const double term_3 = (lambdaT - M*p - c_int) * std::log(arma::det(Sigma_eps)) / 2.0;

    return cons_val + term_1 + term_2 + term_3;
}

template<typename T>
double
dsgevar<T>::log_likelihood_inf(const arma::mat& Gamma_YY, const arma::mat& Gamma_YX, const arma::mat& Gamma_XX)
const
{
    const arma::mat beta = arma::inv_sympd(Gamma_XX)*Gamma_YX.t();
    const arma::mat Sigma_eps = Gamma_YY - Gamma_YX*beta;

    const arma::mat Sigma_bar_eps = YY + beta.t()*XX*beta - beta.t()*XY - XY.t()*beta;

    const double cons_val = - M*(n-p)*std::log(2*arma::datum::pi)/2.0;
    const double logLikelihood = cons_val - (n-p)*( std::log(arma::det(Sigma_eps)) + arma::trace(arma::inv_sympd(Sigma_eps)*Sigma_bar_eps) ) / 2.0;

    return logLikelihood;
}

//
// log posterior kernel

template<typename T>
double
dsgevar<T>::log_posterior_kernel(const arma::vec& pars_inp)
const
{
    const bool finite_lambda = std::isfinite(lambda);

    // setup DSGE model and solve

    dsge_obj.solve_to_state_space(pars_inp);

    // compute model-implied moments

    arma::mat Gamma_YY, Gamma_YX, Gamma_XX;

    model_moments(Gamma_YY, Gamma_YX, Gamma_XX);

    // calculate log-likelihood

    double log_likelihood_val = 0.0;

    if (finite_lambda) {

        const double tau = lambda / (1.0 + lambda);

        const arma::mat Gamma_bar_YY = tau*Gamma_YY + (1.0 - tau)*YY;
        const arma::mat Gamma_bar_YX = tau*Gamma_YX.t() + (1.0 - tau)*XY;
        const arma::mat Gamma_bar_XX = tau*Gamma_XX + (1.0 - tau)*XX;

        log_likelihood_val = log_likelihood(Gamma_YY,Gamma_YX,Gamma_XX,Gamma_bar_YY,Gamma_bar_YX,Gamma_bar_XX);
    } else {
        log_likelihood_val = log_likelihood_inf(Gamma_YY,Gamma_YX,Gamma_XX);
    }

    // compute the prior

    const double log_prior_val = dsge_obj.log_prior(pars_inp);

    //

    return log_likelihood_val + log_prior_val;
}

//
// estimate posterior mode

template<typename T>
double
dsgevar<T>::mode_objfn(const arma::vec& pars_inp, arma::vec* grad_vec, void* mode_data)
{
    dsgevar_estim_data<T>* dta = reinterpret_cast<dsgevar_estim_data<T>*>(mode_data);
    dsgevar<T> dsgevar_obj = dta->dsgevar_obj;

    return - dsgevar_obj.log_posterior_kernel(pars_inp);
}

template<typename T>
arma::vec
dsgevar<T>::estim_mode(const arma::vec& initial_vals)
{
    return this->estim_mode(initial_vals,nullptr,nullptr);
}

template<typename T>
arma::vec
dsgevar<T>::estim_mode(const arma::vec& initial_vals, arma::mat& vcov_mat)
{
    return this->estim_mode(initial_vals,&vcov_mat,nullptr);
}

template<typename T>
arma::vec
dsgevar<T>::estim_mode(const arma::vec& initial_vals, arma::mat* vcov_mat, optim::algo_settings* settings_inp)
{
    dsgevar_estim_data<T> mode_data;
    mode_data.dsgevar_obj = *this;

    //
    // optimization settings

    optim::algo_settings settings;

    if (settings_inp) {
        settings = *settings_inp;
    }

    bool vals_bound = (dsge_obj.lower_bounds.n_elem > 0 || dsge_obj.upper_bounds.n_elem > 0) ? true : false;
    settings.vals_bound = vals_bound;

    settings.lower_bounds = dsge_obj.lower_bounds;
    settings.upper_bounds = dsge_obj.upper_bounds;

    //
    // run optim

    arma::vec ret_vec = initial_vals;

    optim::de(ret_vec,mode_objfn,&mode_data,settings);

    // compute standard errors

    if (vcov_mat) {
        arma::mat hess_mat = optim::numerical_hessian(ret_vec,nullptr,mode_objfn,&mode_data);

        *vcov_mat = arma::inv(hess_mat);
    }

    //

    return ret_vec;
}

//
// MCMC

template<typename T>
double
dsgevar<T>::mcmc_objfn(const arma::vec& pars_inp, void* mcmc_data)
{
    dsgevar_estim_data<T>* dta = reinterpret_cast<dsgevar_estim_data<T>*>(mcmc_data);
    dsgevar<T> dsgevar_obj = dta->dsgevar_obj; // need to create a copy to be thread safe

    return dsgevar_obj.log_posterior_kernel(pars_inp);
}

template<typename T>
void
dsgevar<T>::estim_mcmc(const arma::vec& initial_vals, mcmc::algo_settings_t* settings_inp)
{

    dsgevar_estim_data<T> mcmc_data;
    mcmc_data.dsgevar_obj = *this;

    //
    // MCMC settings

    mcmc::algo_settings_t settings;

    if (settings_inp) {
        settings = *settings_inp;
    }

    bool vals_bound = (dsge_obj.lower_bounds.n_elem > 0 || dsge_obj.upper_bounds.n_elem > 0) ? true : false;
    settings.vals_bound = vals_bound;

    settings.lower_bounds = dsge_obj.lower_bounds;
    settings.upper_bounds = dsge_obj.upper_bounds;

    //
    // run MCMC

    arma::cube res_cube;

    mcmc::de(initial_vals,res_cube,mcmc_objfn,&mcmc_data,settings);

    dsge_obj.dsge_draws = std::move(cube_to_mat(res_cube));

    gibbs();
}

template<typename T>
void
dsgevar<T>::gibbs()
{
    const bool finite_lambda = std::isfinite(lambda);

    const int n_draws = dsge_obj.dsge_draws.n_rows;

    beta_draws.set_size(K, M, n_draws);
    Sigma_draws.set_size(M, M, n_draws);

    Sigma_pt_dof = (1+lambda)*(n-p) - M*p - c_int; // // will be inf (?) in the case of lambda = Inf

    //
    // begin loop

    dsge<T> dsge_obj_copy = dsge_obj; // thread safety

    dsge_obj_copy.estim_data.reset();
    dsge_obj_copy.dsge_draws.reset();

#ifdef BM_USE_OPENMP
    #pragma omp parallel for firstprivate(dsge_obj_copy)
#endif
    for (int i=0; i < n_draws; i++) {

        // access saved DSGE parameters and solve the model

        dsge_obj_copy.solve_to_state_space(dsge_obj.dsge_draws.row(i).t());

        // get the model-implied moments

        arma::mat Gamma_YY, Gamma_YX, Gamma_XX;

        model_moments(dsge_obj_copy, Gamma_YY, Gamma_YX, Gamma_XX);

        arma::mat beta_draw, Sigma_draw;

        if (finite_lambda) {

            const double tau = lambda / (1.0 + lambda); // will be NaN in the case of lambda = Inf

            arma::mat Gamma_bar_YY = tau*Gamma_YY + (1.0 - tau)*YY;
            arma::mat Gamma_bar_YX = tau*Gamma_YX.t() + (1.0 - tau)*XY;
            arma::mat Gamma_bar_XX = tau*Gamma_XX + (1.0 - tau)*XX;

            // draw Sigma

            arma::mat beta_hat = arma::inv_sympd(Gamma_bar_XX)*Gamma_bar_YX;

            arma::mat Sigma_eps = Gamma_bar_YY - Gamma_bar_YX.t()*beta_hat;

            Sigma_draw = stats::rinvwish<arma::mat,double,int>((1+lambda)*(n-p)*Sigma_eps, Sigma_pt_dof);

            // draw beta

            arma::mat alpha_pt_var_b  = arma::kron(Sigma_eps,arma::inv_sympd(lambda*Gamma_XX + XX) / (double) (n-p));

            beta_draw = arma::reshape( stats::rmvnorm<arma::mat>(arma::vectorise(beta_hat), alpha_pt_var_b), K,M);

        } else {

            beta_draw = arma::inv_sympd(Gamma_XX)*Gamma_YX.t();
            Sigma_draw = Gamma_YY - Gamma_YX*beta_draw;

        }

        // save draws

        beta_draws.slice(i)  = std::move(beta_draw);
        Sigma_draws.slice(i) = std::move(Sigma_draw);
    }
    //
    alpha_pt_mean = arma::vectorise(arma::mean(beta_draws,2));
    Sigma_pt_mean = arma::mean(Sigma_draws,2);

    alpha_pt_var = arma::cov( cube_to_mat(beta_draws,true) );
    //
}

//
// IRFs

template<typename T>
arma::cube
dsgevar<T>::IRF(const int n_irf_periods)
const
{
    const int n_draws = beta_draws.n_slices;
    // const int K_adj = K - n_ext_vars;
    const int K_adj = K;

    arma::cube irfs_ret(M, M, n_irf_periods*n_draws);

    //

    dsge<T> dsge_obj_copy = dsge_obj; // thread safety

    dsge_obj_copy.estim_data.reset();
    dsge_obj_copy.dsge_draws.reset();

#ifdef BM_USE_OPENMP
    #pragma omp parallel for firstprivate(dsge_obj_copy)
#endif
    for (int j=1; j <= n_draws; j++) {
        arma::mat beta_b = beta_draws(arma::span(c_int,K_adj-1),arma::span(),arma::span(j-1,j-1)); // b'th draw, minus coefficients on any external variables

        arma::mat Sigma_b = Sigma_draws.slice(j-1);
        arma::mat impact_mat = arma::chol(Sigma_b,"lower");

        // adjust the impact matrix

        arma::mat shocks_cov, G_state;
        dsge_obj_copy.model_fn(dsge_obj.dsge_draws.row(j-1).t(),dsge_obj_copy.lrem_obj,shocks_cov,dsge_obj_copy.kalman_mat_C,dsge_obj_copy.kalman_mat_H,dsge_obj_copy.kalman_mat_R);

        dsge_obj_copy.solve();

        dsge_obj_copy.state_space(dsge_obj_copy.kalman_mat_F,G_state);

        arma::mat dsge_obs_impact = dsge_obj_copy.kalman_mat_H.t() * G_state*arma::chol(shocks_cov,"lower");

        arma::mat Q_dcm, R_dcm;
        arma::qr(Q_dcm,R_dcm,dsge_obs_impact.t());

        Q_dcm *= -1.0;
        R_dcm *= -1.0;

        arma::vec R_diag = R_dcm.diag();
        arma::uvec R_pos_inds = arma::find(R_diag > 0.0);

        arma::vec S = - arma::ones((int) dsge_obj_copy.kalman_mat_H.n_cols,1);

        if (R_pos_inds.n_elem > 0) {
            S.rows(R_pos_inds).fill( 1.0 );
        }

        arma::mat Q = Q_dcm*arma::diagmat(S);

        impact_mat = impact_mat*Q.t();

        //

        arma::mat impact_mat_b = arma::zeros(K_adj-c_int,M);
        arma::mat impact_mat_h = arma::zeros(M,M);

        impact_mat_b.rows(0,M-1) = impact_mat;

        irfs_ret.slice((j-1)*n_irf_periods) = impact_mat;

        for (int i=2; i <= n_irf_periods; i++) {
            impact_mat_h = beta_b.t()*impact_mat_b;
            irfs_ret.slice((j-1)*n_irf_periods + (i-1)) = impact_mat_h;

            if(K_adj > M + c_int){
                impact_mat_b.rows(M,K_adj-c_int-1) = impact_mat_b.rows(0,K_adj-M-c_int-1);
            }

            impact_mat_b.rows(0,M-1) = impact_mat_h;
        }
    }

    //

    return irfs_ret;
}

//
// forecasting

template<typename T>
arma::cube
dsgevar<T>::forecast(const int horizon, const bool incl_shocks)
{
    return this->forecast_int(nullptr,horizon,incl_shocks);
}

template<typename T>
arma::cube
dsgevar<T>::forecast(const arma::mat& X_T, const int horizon, const bool incl_shocks)
{
    return this->forecast_int(&X_T,horizon,incl_shocks);
}

template<typename T>
arma::cube
dsgevar<T>::forecast_int(const arma::mat* X_T_inp, const int horizon, const bool incl_shocks)
{
    const int n_draws = beta_draws.n_slices;
    const int K_adj = K;

    arma::mat beta_b(K_adj,M), Sigma_b(M,M);       // bth draw

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
    
    arma::mat X_Th = X_T;

    arma::mat Y_forecast(horizon,M);
    arma::cube forecast_mat(horizon, M, n_draws);
    //
    if (incl_shocks) {
        arma::mat chol_shock_cov;

        for (int i=0; i < n_draws; i++) {
            beta_b  = beta_draws.slice(i);
            Sigma_b = Sigma_draws.slice(i);

            chol_shock_cov = arma::chol(Sigma_b,"lower");

            Y_forecast.zeros();
            X_Th = X_T;

            for (int j=1; j<=horizon; j++) {
                Y_forecast.row(j-1) = X_Th*beta_b + arma::trans(stats::rmvnorm<arma::mat>(arma::zeros(Sigma_b.n_rows,1),chol_shock_cov,true));

                if (K_adj > M + c_int) {
                    X_Th(0,arma::span(M+c_int,K_adj-1)) = X_Th(0,arma::span(c_int,K_adj-M-1));
                }

                X_Th(0,arma::span(c_int,M-1+c_int)) = Y_forecast.row(j-1);
            }
            //
            forecast_mat.slice(i) = Y_forecast;
        }
    } else {
        for (int i=0; i < n_draws; i++) {
            beta_b = beta_draws.slice(i);

            Y_forecast.zeros();
            X_Th = X_T;

            for (int j=1; j <= horizon; j++) {
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

//
// get filtered states

template<typename T>
arma::cube
dsgevar<T>::state_filter()
const
{
    return dsge_obj.state_filter();
}
