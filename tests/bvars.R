#
# BVAR with steady-state prior
# 
rm(list=ls())
library(BMR)
#
data(BMRVARData)
USMacroData<-USMacroData[,2:4]
#
mypsi <- c(3,6,5)
mycfp <- c(0.9,0.9,0.9)
#
# Different p
#
testbvars <- BVARS(USMacroData,mypsi,mycfp,p=1,irf.periods=20,keep=20000,burnin=5000,XiPsi=1,HP1=0.5,HP4=2,gamma=NULL)
plot(testbvars,save=F)
IRF(testbvars,save=F)
forecast(testbvars,shocks=T,backdata=10,save=F)
#
testbvars <- BVARS(USMacroData,mypsi,mycfp,p=2,irf.periods=20,keep=20000,burnin=5000,XiPsi=1,HP1=0.5,HP4=2,gamma=NULL)
plot(testbvars,save=F)
IRF(testbvars,save=F)
forecast(testbvars,shocks=T,backdata=10,save=F)
#
testbvars <- BVARS(USMacroData,mypsi,mycfp,p=3,irf.periods=20,keep=20000,burnin=5000,XiPsi=1,HP1=0.5,HP4=2,gamma=NULL)
plot(testbvars,save=F)
IRF(testbvars,save=F)
forecast(testbvars,shocks=T,backdata=10,save=F)
#
testbvars <- BVARS(USMacroData,mypsi,mycfp,p=4,irf.periods=20,keep=20000,burnin=5000,XiPsi=1,HP1=0.5,HP4=2,gamma=NULL)
plot(testbvars,save=F)
IRF(testbvars,save=F)
forecast(testbvars,shocks=T,backdata=10,save=F)
#
#
#END