#
# BVAR with time-varying parameters
# 
rm(list=ls())
library(BMR)
#
data(BMRVARData)
USMacroData<-USMacroData[,2:4]
#
irf.points <- c(1979,1988,1997,2006)
yearlab <- seq(1955.00,2010.75,0.25)
#
# Different p
#
bvartvptest <- BVARTVP(USMacroData,yearlab,tau=80,p=1,irf.periods=20,irf.points=irf.points,
                       keep=30000,burnin=40000,XiBeta=4,XiQ=0.005,gammaQ=NULL,XiSigma=1,gammaS=NULL)
plot(bvartvptest,percentiles=c(.16,.50,.84),save=F)
IRF(bvartvptest,percentiles=c(.16,.50,.84),save=F)
#
bvartvptest <- BVARTVP(USMacroData,yearlab,tau=80,p=2,irf.periods=20,irf.points=irf.points,
                       keep=30000,burnin=40000,XiBeta=4,XiQ=0.005,gammaQ=NULL,XiSigma=1,gammaS=NULL)
plot(bvartvptest,percentiles=c(.16,.50,.84),save=F)
IRF(bvartvptest,percentiles=c(.16,.50,.84),save=F)
#
bvartvptest <- BVARTVP(USMacroData,yearlab,tau=80,p=3,irf.periods=20,irf.points=irf.points,
                       keep=30000,burnin=40000,XiBeta=4,XiQ=0.005,gammaQ=NULL,XiSigma=1,gammaS=NULL)
plot(bvartvptest,percentiles=c(.16,.50,.84),save=F)
IRF(bvartvptest,percentiles=c(.16,.50,.84),save=F)
#
bvartvptest <- BVARTVP(USMacroData,yearlab,tau=80,p=4,irf.periods=20,irf.points=irf.points,
                       keep=30000,burnin=40000,XiBeta=4,XiQ=0.005,gammaQ=NULL,XiSigma=1,gammaS=NULL)
plot(bvartvptest,percentiles=c(.16,.50,.84),save=F)
IRF(bvartvptest,percentiles=c(.16,.50,.84),save=F)
#
# 4 variables
#
rm(list=ls())
data(BMRVARData)
USMacroData<-USMacroData[,2:4]
USMacroData<-data.frame(USMacroData,USMacroData[,1]*USMacroData[,1] + rnorm(231,0,0.2))
colnames(USMacroData) <- c(colnames(USMacroData)[1:3],"Keith")
#
irf.points<-c(1979,1988,1997,2006)
yearlab<-seq(1955.00,2010.75,0.25)
#
bvartvptest <- BVARTVP(USMacroData,yearlab,tau=80,p=1,irf.periods=20,irf.points=irf.points,keep=30000,burnin=40000,XiBeta=4,XiQ=0.005,gammaQ=NULL,XiSigma=1,gammaS=NULL)
plot(bvartvptest,percentiles=c(.16,.50,.84),save=F)
IRF(bvartvptest,percentiles=c(.16,.50,.84),save=F)
#
#
#END