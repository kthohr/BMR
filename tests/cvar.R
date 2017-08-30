#
# Estimate a VAR with BMR
# 
rm(list=ls())
library(BMR)
#
data(BMRVARData)
USMacroData<-USMacroData[,2:4]
#
# Different p
#
testcvar <- CVAR(USMacroData,p=1,constant=T,boot=10000)
testcvar <- suppressMessages(CVAR(USMacroData,p=1,constant=T,boot=10000))
IRF(testcvar,save=F)
forecast(testcvar,confint=0.99,backdata=10,save=F)
#
testcvar <- CVAR(USMacroData,p=2,constant=T,boot=10000)
IRF(testcvar,save=F)
forecast(testcvar,confint=0.99,backdata=10,save=F)
#
testcvar <- CVAR(USMacroData,p=3,constant=T,boot=10000)
IRF(testcvar,save=F)
forecast(testcvar,confint=0.99,backdata=10,save=F)
#
testcvar <- CVAR(USMacroData,p=4,constant=T,boot=10000)
IRF(testcvar,save=F)
forecast(testcvar,confint=0.99,backdata=10,save=F)
#
#
#END