#
# Estimate a VAR with BMR
# 

rm(list=ls())
library(BMR.Rcpp)

#

data(BMRVARData)
USMacroData <- USMacroData[,2:4]
var_data <- matrix(c(as.matrix(USMacroData)),ncol=3)

var_obj = new(cvar)

#
# Different p

# p = 1

var_obj$build(var_data,TRUE,1)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(var_obj,varnames=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,varnames=colnames(USMacroData),backdata=10,save=FALSE)

# p = 2

var_obj$reset_draws()
var_obj$build(var_data,TRUE,2)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(var_obj,varnames=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,varnames=colnames(USMacroData),backdata=10,save=FALSE)

# p = 3

var_obj$reset_draws()
var_obj$build(var_data,TRUE,3)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(var_obj,varnames=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,varnames=colnames(USMacroData),backdata=10,save=FALSE)

# p = 4

var_obj$reset_draws()
var_obj$build(var_data,TRUE,4)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(var_obj,varnames=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,varnames=colnames(USMacroData),backdata=10,save=FALSE)

#
#END
