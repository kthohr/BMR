#
# Estimate a VAR with BMR
# 

rm(list=ls())
library(BMR)

#

data(BMRVARData)
var_data <- data.matrix(USMacroData[,2:4])

var_obj <- new(cvar)

#
# Different p

# p = 1

var_obj$build(var_data,TRUE,1)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(var_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)
FEVD(var_obj,20,var_names=colnames(var_data),save=FALSE)

# p = 2

var_obj$reset_draws()
var_obj$build(var_data,TRUE,2)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(var_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)

# p = 3

var_obj$reset_draws()
var_obj$build(var_data,TRUE,3)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(var_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)

# p = 4

var_obj$reset_draws()
var_obj$build(var_data,TRUE,4)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(var_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)

#
#END
