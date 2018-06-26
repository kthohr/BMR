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
 * DSGE model class
 */

//
// settings bounds and priors

template<typename T>
void
dsge<T>::set_bounds(const arma::vec& lower_bounds_inp, const arma::vec& upper_bounds_inp)
{
    lower_bounds = lower_bounds_inp;
    upper_bounds = upper_bounds_inp;
}

template<typename T>
void
dsge<T>::set_prior(const arma::uvec& prior_form_inp, const arma::mat& prior_pars_inp)
{
    prior_form = prior_form_inp;
    prior_pars = prior_pars_inp;
}

//
// wrappers to LREM functions

template<typename T>
void
dsge<T>::solve()
const
{
    lrem_obj.solve();
}

template<typename T>
void
dsge<T>::state_space(arma::mat& F_state, arma::mat& G_state)
const
{
    lrem_obj.state_space(F_state,G_state);
}

template<typename T>
arma::mat
dsge<T>::simulate(const int sim_periods, const int burnin)
const
{
    return lrem_obj.simulate(sim_periods, burnin);
}

template<typename T>
void
dsge<T>::solve_to_state_space(const arma::vec& pars_inp)
const
{
    arma::mat shocks_cov;
    model_fn(pars_inp,lrem_obj,shocks_cov,kalman_mat_C,kalman_mat_H,kalman_mat_R);

    solve();

    arma::mat G_state;
    state_space(kalman_mat_F,G_state);

    kalman_mat_Q = G_state*shocks_cov*G_state.t();

    lrem_obj.shocks_cov = std::move(shocks_cov);
}

//
// prior

template<typename T>
double
dsge<T>::log_prior(const arma::vec& pars_inp)
const
{
    const int n_param = pars_inp.n_elem;

    //

    double log_prior_val = 0;

    for (int i = 0; i < n_param; i++) {
        switch (prior_form(i)) {
            case 1: // normal prior
                log_prior_val += stats::dnorm(pars_inp(i),prior_pars(i,0),prior_pars(i,1),true);
                break;

            case 2: // gamma prior
                log_prior_val += stats::dgamma(pars_inp(i),prior_pars(i,0),prior_pars(i,1),true);
                break;

            case 3: // inverse gamma prior
                log_prior_val += stats::dinvgamma(pars_inp(i),prior_pars(i,0),prior_pars(i,1),true);
                break;

            case 4: // beta prior
                log_prior_val += stats::dbeta(pars_inp(i),prior_pars(i,0),prior_pars(i,1),true);
                break;

            case 5: // uniform prior
                log_prior_val += stats::dunif(pars_inp(i),prior_pars(i,0),prior_pars(i,1),true);
                break;

            default:
                printf("unrecognized prior form.\n");
                break;

        }
    }
    //
    return log_prior_val;
}

//
// posterior kernel

template<typename T>
double
dsge<T>::log_posterior_kernel(const arma::vec& pars_inp)
const
{
    // setup model and solve

    solve_to_state_space(pars_inp);

    // run filter

    double log_likelihood_val = 0.0;

    switch (filter_choice) {
        case 1:
            log_likelihood_val = kalman_filter(estim_data, kalman_mat_F,kalman_mat_Q, kalman_mat_C,kalman_mat_H,kalman_mat_R);
            break;

        case 2:
            log_likelihood_val = chand_recur(estim_data, kalman_mat_F,kalman_mat_Q, kalman_mat_C,kalman_mat_H,kalman_mat_R);
            break;

        default:
            printf("error: unknown choice for filter\n");
            break;

    }

    // compute the prior

    const double log_prior_val = log_prior(pars_inp);

    //
    return log_likelihood_val + log_prior_val;
}

//
// estimate posterior mode

template<typename T>
double
dsge<T>::mode_objfn(const arma::vec& pars_inp, arma::vec* grad_vec, void* mode_data)
{
    dsge_estim_data<T>* dta = reinterpret_cast<dsge_estim_data<T>*>(mode_data);
    dsge<T> dsge_obj = dta->dsge_obj; // need to create a copy to be thread safe

    return - dsge_obj.log_posterior_kernel(pars_inp);
}

template<typename T>
arma::vec
dsge<T>::estim_mode(const arma::vec& initial_vals)
{
    return this->estim_mode(initial_vals,nullptr,nullptr);
}

template<typename T>
arma::vec
dsge<T>::estim_mode(const arma::vec& initial_vals, arma::mat& vcov_mat)
{
    return this->estim_mode(initial_vals,&vcov_mat,nullptr);
}

template<typename T>
arma::vec
dsge<T>::estim_mode(const arma::vec& initial_vals, arma::mat* vcov_mat, optim::algo_settings* settings_inp)
{
    dsge_estim_data<T> mode_data;
    mode_data.dsge_obj = *this;

    //
    // optimization settings

    optim::algo_settings settings;

    if (settings_inp) {
        settings = *settings_inp;
    }

    bool vals_bound = (lower_bounds.n_elem > 0 || upper_bounds.n_elem > 0) ? true : false;
    settings.vals_bound = vals_bound;

    settings.lower_bounds = lower_bounds;
    settings.upper_bounds = upper_bounds;

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
dsge<T>::mcmc_objfn(const arma::vec& pars_inp, void* mcmc_data)
{
    dsge_estim_data<T>* dta = reinterpret_cast<dsge_estim_data<T>*>(mcmc_data);
    dsge<T> dsge_obj = dta->dsge_obj; // need to create a copy to be thread safe

    return dsge_obj.log_posterior_kernel(pars_inp);
}

template<typename T>
void
dsge<T>::estim_mcmc(const arma::vec& initial_vals, mcmc::algo_settings_t* settings_inp)
{

    dsge_estim_data<T> mcmc_data;
    mcmc_data.dsge_obj = *this;

    //
    // MCMC settings

    mcmc::algo_settings_t settings;

    if (settings_inp) {
        settings = *settings_inp;
    }

    bool vals_bound = (lower_bounds.n_elem > 0 || upper_bounds.n_elem > 0) ? true : false;
    settings.vals_bound = vals_bound;

    settings.lower_bounds = lower_bounds;
    settings.upper_bounds = upper_bounds;

    //
    // run MCMC

    arma::cube res_cube;

    mcmc::de(initial_vals,res_cube,mcmc_objfn,&mcmc_data,settings);

    dsge_draws = std::move(cube_to_mat(res_cube));
}

//
// IRFs

template<typename T>
arma::cube
dsge<T>::IRF(const int n_irf_periods, const bool observ_irfs)
const
{
    const int n_draws = dsge_draws.n_rows;

    solve_to_state_space(dsge_draws.row(0).t());

    arma::cube test_cube = lrem_obj.IRF(n_irf_periods);

    const int n_resp_cols = (observ_irfs) ? kalman_mat_H.n_cols : test_cube.n_cols;
    const int n_shocks = test_cube.n_slices;

    arma::cube irfs_ret(test_cube.n_rows, n_resp_cols, n_shocks*n_draws);

    //

    dsge<T> dsge_obj_copy = *this; // thread safety

    dsge_obj_copy.estim_data.reset();
    dsge_obj_copy.dsge_draws.reset();

#ifdef BM_USE_OPENMP
    #pragma omp parallel for firstprivate(dsge_obj_copy)
#endif
    for (int j=0; j < n_draws; j++)
    {    
        dsge_obj_copy.solve_to_state_space(dsge_draws.row(j).t());

        arma::cube irf_mats = dsge_obj_copy.lrem_obj.IRF(n_irf_periods);

        if (observ_irfs) {
            for (int k=0; k < n_shocks; k++) {
                irf_mats.slice(k) = irf_mats.slice(k) * dsge_obj_copy.kalman_mat_H;
            }
        }
        
        irfs_ret.slices(j*n_shocks,(j+1)*n_shocks-1) = irf_mats;

    }

    //

    return irfs_ret;
}

//
// forecasting

template<typename T>
arma::cube
dsge<T>::forecast(const int horizon, const bool incl_shocks)
const
{
    const int M = estim_data.n_cols;
    const int n_draws = dsge_draws.n_rows;

    solve_to_state_space(dsge_draws.row(0).t());

    const int n_states = kalman_mat_F.n_cols;

    arma::cube forecast_cube(horizon,M,n_draws);

    //

    dsge<T> dsge_obj_copy = *this; // thread safety

    dsge_obj_copy.estim_data.reset();
    dsge_obj_copy.dsge_draws.reset();

#ifdef BM_USE_OPENMP
    #pragma omp parallel for firstprivate(dsge_obj_copy)
#endif
    for (int j=0; j < n_draws; j++) {
        
        dsge_obj_copy.solve_to_state_space(dsge_draws.row(j).t());

        arma::mat G_state;
        dsge_obj_copy.state_space(dsge_obj_copy.kalman_mat_F,G_state);

        //

        arma::mat filt_states;

        kalman_filter(estim_data, dsge_obj_copy.kalman_mat_F,dsge_obj_copy.kalman_mat_Q, dsge_obj_copy.kalman_mat_C,dsge_obj_copy.kalman_mat_H,dsge_obj_copy.kalman_mat_R, &filt_states);

        arma::vec state_n = filt_states.row((int) filt_states.n_rows - 1).t();

        //

        arma::mat chol_shocks_cov = arma::chol(dsge_obj_copy.lrem_obj.shocks_cov);

        arma::mat iter_mat = arma::eye(n_states,n_states), forecast_mat(horizon,M);

        for (int i=0; i < horizon; i++) {
            state_n = dsge_obj_copy.kalman_mat_F * state_n;

            if (incl_shocks) {
                state_n += G_state*stats::rmvnorm<arma::mat>(arma::zeros(chol_shocks_cov.n_rows,1),chol_shocks_cov,true);
            }

            forecast_mat.row(i) = arma::trans( dsge_obj_copy.kalman_mat_C + dsge_obj_copy.kalman_mat_H.t() * state_n );
        }

        //
        
        forecast_cube.slice(j) = std::move(forecast_mat);
    }

    //

    return forecast_cube;
}

//
// get filtered states

template<typename T>
arma::cube
dsge<T>::state_filter()
const
{
    const int n = estim_data.n_rows;
    const int n_draws = dsge_draws.n_rows;

    solve_to_state_space(dsge_draws.row(0).t());

    const int n_states = kalman_mat_F.n_cols;

    arma::cube filter_cube(n,n_states,n_draws);

    //

    dsge<T> dsge_obj_copy = *this; // thread safety

    dsge_obj_copy.estim_data.reset();
    dsge_obj_copy.dsge_draws.reset();

#ifdef BM_USE_OPENMP
    #pragma omp parallel for firstprivate(dsge_obj_copy)
#endif
    for (int j=0; j < n_draws; j++) {
        
        dsge_obj_copy.solve_to_state_space(dsge_draws.row(j).t());

        //

        arma::mat filt_vals;

        kalman_filter(estim_data, dsge_obj_copy.kalman_mat_F,dsge_obj_copy.kalman_mat_Q, dsge_obj_copy.kalman_mat_C,dsge_obj_copy.kalman_mat_H,dsge_obj_copy.kalman_mat_R, &filt_vals);

        //

        filter_cube.slice(j) = std::move(filt_vals);
    }

    //

    return filter_cube;
}
