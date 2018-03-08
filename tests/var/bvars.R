#
# BVAR with steady-state prior
# 

rm(list=ls())
library(BMR)

#

data(BMRVARData)
bvar_data <- data.matrix(USMacroData[,2:4])

#

coef_prior <- c(0.9,0.9,0.9)
psi_prior <- matrix(c(3,6,5),nrow=1)

HP1 <- 0.5
HP2 <- 0.5
XiPsi <- 1
gamma = 4

bvar_obj <- new(bvars)

#
# Different p

# p = 1

bvar_obj$build(bvar_data,TRUE,1)
bvar_obj$prior(coef_prior,HP1,HP2,psi_prior,XiPsi,gamma,FALSE)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(bvar_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)
FEVD(bvar_obj,20,var_names=colnames(bvar_data),save=FALSE)

# p = 2

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,2)
bvar_obj$prior(coef_prior,HP1,HP2,psi_prior,XiPsi,gamma,FALSE)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(bvar_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)

# p = 3

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,3)
bvar_obj$prior(coef_prior,HP1,HP2,psi_prior,XiPsi,gamma,FALSE)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(bvar_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)

# p = 4

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,4)
bvar_obj$prior(coef_prior,HP1,HP2,psi_prior,XiPsi,gamma,FALSE)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(bvar_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)

#
#END

