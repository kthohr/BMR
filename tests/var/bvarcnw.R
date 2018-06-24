#
# BVAR with conditional normal-inverse-Wishart prior
# 

rm(list=ls())
library(BMR)

#

data(BMRVARData)
bvar_data <- data.matrix(USMacroData[,2:4])

#

coef_prior <- c(0.9,0.9,0.9)
HP_1 <- 1/2
HP_3 <- 1
gamma <- 4

bvar_obj <- new(bvarcnw)

#
# Different p

# p = 1

bvar_obj$build(bvar_data,TRUE,1)
bvar_obj$prior(coef_prior,HP_1,HP_3,gamma)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,var_names=colnames(bvar_data),save=FALSE)
plot(bvar_obj,var_names=colnames(bvar_data),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(bvar_data),back_data=10,save=FALSE)
FEVD(bvar_obj,20,var_names=colnames(bvar_data),save=FALSE)

# p = 2

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,2)
bvar_obj$prior(coef_prior,HP_1,HP_3,gamma)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,var_names=colnames(bvar_data),save=FALSE)
plot(bvar_obj,var_names=colnames(bvar_data),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(bvar_data),back_data=10,save=FALSE)

# p = 3

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,3)
bvar_obj$prior(coef_prior,HP_1,HP_3,gamma)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,var_names=colnames(bvar_data),save=FALSE)
plot(bvar_obj,var_names=colnames(bvar_data),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(bvar_data),back_data=10,save=FALSE)

# p = 4

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,4)
bvar_obj$prior(coef_prior,HP_1,HP_3,gamma)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,var_names=colnames(bvar_data),save=FALSE)
plot(bvar_obj,var_names=colnames(bvar_data),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(bvar_data),back_data=10,save=FALSE)

#
#END
