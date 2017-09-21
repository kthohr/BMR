#
# BVAR with time-varying parameters

rm(list=ls())
library(BMR)

#

data(BMRVARData)
bvar_data <- data.matrix(USMacroData[,2:4])

#

tau <- 80

XiBeta <- 4
XiQ <- 0.001
gammaQ <- tau
XiSigma <- 1
gammaS = 4

which_irfs = c(17,53,89,125)

bvar_obj <- new(bvartvp)

#
# Different p

# p = 1

bvar_obj$build(bvar_data,TRUE,1)
bvar_obj$prior(tau,XiBeta,XiQ,gammaQ,XiSigma,gammaS)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,which_irfs,var_names=colnames(USMacroData))
plot(bvar_obj,var_names=colnames(USMacroData))

# p = 2

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,2)
bvar_obj$prior(tau,XiBeta,XiQ,gammaQ,XiSigma,gammaS)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,which_irfs,var_names=colnames(USMacroData))
plot(bvar_obj,var_names=colnames(USMacroData))

# p = 3

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,3)
bvar_obj$prior(tau,XiBeta,XiQ,gammaQ,XiSigma,gammaS)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,which_irfs,var_names=colnames(USMacroData))
plot(bvar_obj,var_names=colnames(USMacroData))

# p = 4

bvar_obj$reset_draws()
bvar_obj$build(bvar_data,TRUE,4)
bvar_obj$prior(tau,XiBeta,XiQ,gammaQ,XiSigma,gammaS)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,which_irfs,var_names=colnames(USMacroData))
plot(bvar_obj,var_names=colnames(USMacroData))

#
#END
