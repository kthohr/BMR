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
  ##   You should have received a copy of the GNU General Public License
  ##   along with BMLib. If not, see <http://www.gnu.org/licenses/>.
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
{
    lrem_obj.solve();
}

template<typename T>
void
dsge<T>::state_space(arma::mat& F_state, arma::mat& G_state)
{
    lrem_obj.state_space(F_state,G_state);
}

template<typename T>
arma::mat
dsge<T>::simulate(const int sim_periods, const int burnin)
{
    return lrem_obj.simulate(sim_periods, burnin);
}

template<typename T>
void
dsge<T>::solve_to_state_space(const arma::vec& pars_inp)
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
{
    // setup model and solve

    solve_to_state_space(pars_inp);

    // run filter

    const double log_likelihood_val = kalman_filter(estim_data, kalman_mat_F,kalman_mat_Q, kalman_mat_C,kalman_mat_H,kalman_mat_R);

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
    return this->estim_mode(initial_vals,nullptr);
}

template<typename T>
arma::vec
dsge<T>::estim_mode(const arma::vec& initial_vals, optim::opt_settings* settings_inp)
{
    dsge_estim_data<T> mode_data;
    mode_data.dsge_obj = *this;

    //
    // optimization settings

    optim::opt_settings settings;

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
dsge<T>::estim_mcmc(const arma::vec& initial_vals, mcmc::mcmc_settings* settings_inp)
{

    dsge_estim_data<T> mcmc_data;
    mcmc_data.dsge_obj = *this;

    //
    // MCMC settings

    mcmc::mcmc_settings settings;

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
dsge<T>::IRF(const int n_irf_periods)
{
    const int n_draws = dsge_draws.n_rows;

    solve_to_state_space(dsge_draws.row(0).t());

    arma::cube test_cube = lrem_obj.IRF(n_irf_periods);

    const int n_shocks = test_cube.n_slices;

    arma::cube irfs(test_cube.n_rows, test_cube.n_cols, n_shocks*n_draws);

    //

    dsge<T> dsge_obj_copy = *this; // thread safety

    dsge_obj_copy.estim_data.reset();
    dsge_obj_copy.dsge_draws.reset();

#ifdef BM_USE_OMP
    #pragma omp parallel for firstprivate(dsge_obj_copy)
#endif
    for (int j=0; j < n_draws; j++) {
        
        dsge_obj_copy.solve_to_state_space(dsge_draws.row(j).t());

        irfs.slices(j*n_shocks,(j+1)*n_shocks-1) = dsge_obj_copy.lrem_obj.IRF(n_irf_periods);

    }

    //

    return irfs;
}
