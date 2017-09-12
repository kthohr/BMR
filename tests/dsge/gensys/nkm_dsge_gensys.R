
#

rm(list=ls())

library(BMR)

source("nkm_model.R")

#

data(BMRVARData)
dsgedata <- USMacroData[24:211,-c(1,3)]
dsgedata <- as.matrix(dsgedata)
for(i in 1:2){
    dsgedata[,i] <- dsgedata[,i] - mean(dsgedata[,i])
}

#

obj <- new(dsge_gensys)
obj$set_model_fn(nkm_model_simple)

x <- c(1)

obj$eval_model(x)

lrem_obj = obj$lrem
lrem_obj$solve()

lrem_obj$shocks_cov <- matrix(c(1,0,0,0.125),2,2,byrow=TRUE)

sim_data <- lrem_obj$simulate(200,800)$sim_vals

sim_data <- cbind(sim_data[,3],sim_data[,5])

#

prior_pars <- cbind(c(1.0),
                    c(0.05))

prior_form <- c(1)

obj$set_prior(prior_form,prior_pars)

#

par_bounds <- cbind(c(-Inf),
                    c( Inf))

opt_bounds <- cbind(c(0.5),
                    c(3.0))

obj$set_bounds(opt_bounds[,1],opt_bounds[,2])

obj$opt_initial_lb <- opt_bounds[,1]
obj$opt_initial_ub <- opt_bounds[,2]

#

obj$estim_data = sim_data;

#obj$estim_mode(x)

obj$mcmc_initial_lb <- opt_bounds[,1]
obj$mcmc_initial_ub <- opt_bounds[,2]

obj$estim_mcmc(x,50,100,100)

varnames <- c("Output Gap","Output","Inflation","Natural Int","Nominal Int","Labour Supply",
              "Technology","MonetaryPolicy")

plot(obj,parnames="eta",save=FALSE)
IRF(obj,20,varnames=varnames,save=FALSE)
