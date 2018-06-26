# BMR &nbsp; [![Build Status](https://travis-ci.org/kthohr/BMR.svg)](https://travis-ci.org/kthohr/BMR) [![Build status](https://ci.appveyor.com/api/projects/status/github/kthohr/BMR?branch=master)](https://ci.appveyor.com/project/kthohr/BMR/branch/master)


Bayesian Macroeconometrics in R (BMR) is an R interface to BM++, a templated C++ library for estimating Bayesian Vector Autoregression (BVAR) and Dynamic Stochastic General Equilibrium (DSGE) models.

Features:

* BVAR models with flexible prior specification, including: the Minnesota prior; normal-inverse-Wishart prior; and Matias Villani's steady-state prior.
* BVARs with time-varying parameters.
* Solve DSGE models using the method of undetermined coefficients (Uhlig's method) and Sims' QZ method.
* Estimate DSGE and DSGE-VAR models.

## Installation

The simplest installation method is via the devtools package.

```R
install.packages("devtools")
library(devtools)
install_github("kthohr/BMR")
```

Note that BMR requires compilation, so approproate development tools are necessary to install the package.
* Windows users should get [Rtools](https://cran.r-project.org/bin/windows/Rtools/).
* Mac uses should check [here](https://cran.r-project.org/bin/macosx/tools/).

## Documentation

See [http://www.kthohr.com/bmr](http://www.kthohr.com/bmr) for documentation and replication files.

## BVAR Example

Estimate a BVAR model with Minnesota prior:

```R
library("BMR")

data(BMRMCData)

# create bvarm object
bvar_obj <- new(bvarm) 

# add data; add a constant; and set lags (p) = 2
bvar_obj$build(bvarMCdata,TRUE,2)

# prior of 0.5 on own-first-lags, and take default values for the hyperparameters
bvar_obj$prior(c(0.5,0.5))

# draw from the posterior distribution
bvar_obj$gibbs(10000)

# posterior mean and variance
bvar_obj$alpha_pt_mean
bvar_obj$alpha_pt_var

# plot IRFs
IRF(bvar_obj,20,save=FALSE)
```

## DSGE Example

Solve a three-equation New-Keynesian model:

```R
#
# Solve the NK model with gensys
#

rm(list=ls())
library(BMR)

#

alpha    <- 0.33
vartheta <- 6
beta     <- 0.99
theta    <- 0.6667

eta    <- 1               
phi    <- 1                  
phi_pi <- 1.5             
phi_y  <- 0.5/4
rho_a  <- 0.90
rho_v  <- 0.5

BigTheta <- (1-alpha)/(1-alpha+alpha*vartheta)
kappa    <- (((1-theta)*(1-beta*theta))/(theta))*BigTheta*((1/eta)+((phi+alpha)/(1-alpha)))
psi      <- (eta*(1+phi))/(1-alpha+eta*(phi + alpha))

sigma_T <- 1
sigma_M <- 0.25

G0_47 <- (1/eta)*psi*(rho_a - 1)

#Order:                 yg,       y,      pi,      rn,       i,       n,       a,       v,  yg_t+1,  pi_t+1
Gamma0 <- rbind(c(      -1,       0,       0,     eta,  -eta/4,       0,       0,       0,       1,   eta/4),
                c(   kappa,       0,    -1/4,       0,       0,       0,       0,       0,       0,  beta/4),
                c(   phi_y,       0,phi_pi/4,       0,    -1/4,       0,       0,       1,       0,       0),
                c(       0,       0,       0,      -1,       0,       0,   G0_47,       0,       0,       0),
                c(       0,      -1,       0,       0,       0, 1-alpha,       1,       0,       0,       0),
                c(      -1,       1,       0,       0,       0,       0,    -psi,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       1,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       1,       0,       0),
                c(       1,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       1,       0,       0,       0,       0,       0,       0,       0))

Gamma1 <- rbind(c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,   rho_a,       0,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,   rho_v,       0,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       1,       0),
                c(       0,       0,       0,       0,       0,       0,       0,       0,       0,       1))

#

C <- matrix(0,10,1)

Psi <- matrix(0,10,2)
Psi[7,1] <- 1
Psi[8,2] <- 1

Pi <- matrix(0,10,2)
Pi[9,1] <- 1
Pi[10,2] <- 1

#

Sigma <- rbind(c(sigma_T^2,         0),  
               c(        0, sigma_M^2))

#
# build the model and solve it

dsge_obj <- new(gensys)

dsge_obj$build(Gamma0,Gamma1,C,Psi,Pi)
dsge_obj$solve()

# Y_t = G *  Y_{t-1} + Impact * e_t

dsge_obj$G_sol
dsge_obj$impact_sol

# set shocks (Sigma); then plot IRFs and simulate some data

dsge_obj$shocks_cov = Sigma

var_names <- c("Output Gap","Output","Inflation","Natural Int",
               "Nominal Int","Labour Supply","Technology","MonetaryPolicy")

dsge_irf <- IRF(dsge_obj,12,var_names=var_names)
dsge_sim <- dsge_obj$simulate(800,200)
```


## Author

Keith O'Hara

## License

GPL (>= 2) 