/*################################################################################
  ##
  ##   Copyright (C) 2011-2017 Keith O'Hara
  ##
  ##   This file is part of the R package BMR.
  ##
  ##   The R package BMR is free software: you can redistribute it and/or modify
  ##   it under the terms of the GNU General Public License as published by
  ##   the Free Software Foundation, either version 2 of the License, or
  ##   (at your option) any later version.
  ##
  ##   The R package BMR is distributed in the hope that it will be useful,
  ##   but WITHOUT ANY WARRANTY; without even the implied warranty of
  ##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ##   GNU General Public License for more details.
  ##
  ##   You should have received a copy of the GNU General Public License
  ##   along with BMR. If not, see <http://www.gnu.org/licenses/>.
  ##
  ################################################################################*/

class gensys_R : public bm::gensys
{
    public:
        void build_R(arma::mat Gamma_0_inp, arma::mat Gamma_1_inp, arma::mat Gamma_C_inp, arma::mat Psi_inp, arma::mat Pi_inp);
        void solve_R();
        void state_space_R();
        SEXP simulate_R(int n_sim_periods, int n_burnin);
        SEXP IRF_R(int n_irf_periods);
};

class uhlig_R : public bm::uhlig
{
    public:
        void build_R(arma::mat A_inp, arma::mat B_inp, arma::mat C_inp, arma::mat D_inp, arma::mat F_inp, arma::mat G_inp, arma::mat H_inp, arma::mat J_inp, arma::mat K_inp, arma::mat L_inp, arma::mat M_inp, arma::mat N_inp);
        void build_pre_R(arma::mat A_inp, arma::mat B_inp, arma::mat C_inp, arma::mat D_inp);
        void build_exp_R(arma::mat F_inp, arma::mat G_inp, arma::mat H_inp, arma::mat J_inp, arma::mat K_inp, arma::mat L_inp, arma::mat M_inp);
        void build_exog_R(arma::mat N_inp);
        void solve_R();
        void state_space_R();
        SEXP simulate_R(int n_sim_periods, int n_burnin);
        SEXP IRF_R(int n_irf_periods);
};
