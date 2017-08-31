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

#
# Different p

bvar_obj = new(bvarm_R)
bvar_obj$build(bvar_data,TRUE,4)
bvar_obj$prior(coef_prior,1,1,0.5,0.5,100.0,1.0)
bvar_obj$gibbs(10000)
bvar_obj$IRF(20)

IRF(bvar_obj,20,varnames=colnames(USMacroData),save=FALSE)
plot(bvar_obj,varnames=colnames(USMacroData),save=FALSE)


testbvarm <- BVARM(USMacroData,prior,p=1,constant=T,keep=10000,burnin=5000,VType=1,HP1=0.5,HP2=0.5,HP3=100)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=T,backdata=10,save=F)
#
testbvarm <- BVARM(USMacroData,prior,p=2,constant=T,keep=10000,burnin=5000,VType=1,HP1=0.5,HP2=0.5,HP3=100)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=T,backdata=10,save=F)
#
testbvarm <- BVARM(USMacroData,prior,p=3,constant=T,keep=10000,burnin=5000,VType=1,HP1=0.5,HP2=0.5,HP3=100)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=T,backdata=10,save=F)
#
testbvarm <- BVARM(USMacroData,prior,p=4,constant=T,keep=10000,burnin=5000,VType=1,HP1=0.5,HP2=0.5,HP3=100)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=T,backdata=10,save=F)
#
#
#END