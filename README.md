# BMR &nbsp; [![Build Status](https://travis-ci.org/kthohr/BMR.svg)](https://travis-ci.org/kthohr/BMR) [![Build status](https://ci.appveyor.com/api/projects/status/github/kthohr/BMR?branch=master)](https://ci.appveyor.com/project/kthohr/BMR/branch/master)


Bayesian Macroeconometrics in R (BMR) is an R interface to BMLib, a templated C++ library for estimating Bayesian Vector Autoregression (BVAR) and Dynamic Stochastic General Equilibrium (DSGE) models.

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

## Example

Estimate a BVAR model with Minnesota prior.

```R
library("BMR")

data(BMRMCData)

bvar_obj <- new(bvarm) 

bvar_obj$build(bvarMCdata,TRUE,2)
bvar_obj$prior(c(0.5,0.5),1,1,0.5,0.5,100.0,1.0)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,save=FALSE)
```

## Author

Keith O'Hara

## License

GPL (>= 2) 