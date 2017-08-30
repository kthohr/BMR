#
# BVAR with normal-inverse-Wishart prior
# 
rm(list=ls())
library(BMR)
#
data(BMRVARData)
USMacroData<-USMacroData[,2:4]
#
mycfp <- c(0.9,0.9,0.9)
#
# Different p
#
testbvarw <- BVARW(USMacroData,1,mycfp,p=1,constant=T,irf.periods=20,keep=10000,burnin=5000,XiBeta=4,XiSigma=1,gamma=NULL)
testbvarw <- suppressMessages(BVARW(USMacroData,1,mycfp,p=1,constant=T,irf.periods=20,keep=10000,burnin=5000,XiBeta=4,XiSigma=1,gamma=NULL))
plot(testbvarw,save=F)
IRF(testbvarw,save=F)
forecast(testbvarw,periods=10,shocks=T,backdata=10,save=F)
#
testbvarw <- BVARW(USMacroData,1,mycfp,p=2,constant=T,irf.periods=20,keep=10000,burnin=5000,XiBeta=4,XiSigma=1,gamma=NULL)
plot(testbvarw,save=F)
IRF(testbvarw,save=F)
forecast(testbvarw,periods=10,shocks=T,backdata=10,save=F)
#
testbvarw <- BVARW(USMacroData,1,mycfp,p=3,constant=T,irf.periods=20,keep=10000,burnin=5000,XiBeta=4,XiSigma=1,gamma=NULL)
plot(testbvarw,save=F)
IRF(testbvarw,save=F)
forecast(testbvarw,periods=10,shocks=T,backdata=10,save=F)
#
testbvarw <- BVARW(USMacroData,1,mycfp,p=4,constant=T,irf.periods=20,keep=10000,burnin=5000,XiBeta=4,XiSigma=1,gamma=NULL)
plot(testbvarw,save=F)
IRF(testbvarw,save=F)
forecast(testbvarw,periods=10,shocks=T,backdata=10,save=F)
#
# 2 cores
#
testbvarw <- BVARW(USMacroData,cores=2,mycfp,p=4,constant=T,irf.periods=20,keep=10000,burnin=5000,XiBeta=4,XiSigma=1,gamma=NULL)
plot(testbvarw,save=F)
IRF(testbvarw,save=F)
forecast(testbvarw,periods=10,shocks=T,backdata=10,save=F)
#
#
#END