#
# BVAR with Minnesota prior
#

rm(list=ls())
library(BMR.Rcpp)

#

data(BMRVARData)
USMacroData <- USMacroData[,2:4]
bvar_data <- matrix(c(as.matrix(USMacroData)),ncol=3)

#

coef_prior <- c(0.9,0.9,0.9)

bvar_obj = new(bvarm)

#
# Different p

# p = 1

bvar_obj$build(bvar_data,TRUE,1)
bvar_obj$prior(coef_prior,1,1,0.5,0.5,100.0,1.0)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(bvar_obj,varnames=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,varnames=colnames(USMacroData),backdata=10,save=FALSE)

# p = 2

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,2)
bvar_obj$prior(coef_prior,1,1,0.5,0.5,100.0,1.0)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(bvar_obj,varnames=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,varnames=colnames(USMacroData),backdata=10,save=FALSE)

# p = 3

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,3)
bvar_obj$prior(coef_prior,1,1,0.5,0.5,100.0,1.0)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(bvar_obj,varnames=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,varnames=colnames(USMacroData),backdata=10,save=FALSE)

# p = 4

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,4)
bvar_obj$prior(coef_prior,1,1,0.5,0.5,100.0,1.0)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(bvar_obj,varnames=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,varnames=colnames(USMacroData),backdata=10,save=FALSE)

#
#END
