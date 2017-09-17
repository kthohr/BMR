\name{plot}

\alias{plot.bvarm}
\alias{plot.bvars}
\alias{plot.bvartvp}
\alias{plot.bvarw}
\alias{plot.DSGEVAR}
\alias{plot.EDSGE}

\title{Plot BMR Objects.}

\description{
  Plot the output of a BVAR, DSGEVAR, or EDSGE object.
}

\usage{
  \method{plot}{bvarm}(x,type=1,save=FALSE,height=13,width=13,...)
  \method{plot}{bvars}(x,type=1,plotSigma=TRUE,save=FALSE,height=13,width=13,...)
  \method{plot}{bvartvp}(x,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13,...)
  \method{plot}{bvarw}(x,type=1,plotSigma=TRUE,save=FALSE,height=13,width=13,...)
  
  \method{plot}{DSGEVAR}(x,par_names=NULL,BinDenom=40,MCMCplot=FALSE,save=FALSE,height=13,width=13,...)
  \method{plot}{EDSGE}(x,par_names=NULL,BinDenom=40,MCMCplot=FALSE,save=FALSE,height=13,width=13,...)
}

\arguments{
  For objects of class `bvarm',
  \item{x}{
  An object of the above class.}
	\item{type}{
	An integer value indicating the plot style; type=1 will produce a histogram, while type=2 will use smoothed densities.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=FALSE, use this to set the height of the plot.}
	\item{width}{
	If save=FALSE, use this to set the width of the plot.}\cr
  
  For objects of class `bvars' or `bvarw',
  \item{plotSigma}{
  Whether to plot the elements of the residual covariance matrix.}\cr
  
  For objects of class `bvartvp',
  \item{percentiles}{
  Which percentiles of the posterior distribution the function should use.}
  
  For objects of class `DSGEVAR' or `EDSGE',
  \item{par_names}{
  A vector of expressions containing the name of each parameter.}
  \item{BinDenom}{
  Bin width.}
  \item{MCMCplot}{
  Whether to plot a trace of the MCMC run.}\cr
  
  \item{\dots}{
  Additional arguments (not used).}
}

\details{
  For examples of the plotting features, see the accompanying vignette.
}

\value{
  Plots the marginal posterior distributions of the relevant coefficients/parameters.
}

\author{Keith O'Hara}

\examples{
\dontrun{
data(BMRVARData)
prior<-c(0.9,0.95,0.95)
testbvarm <- bvarm(USMacroData[,2:4],prior,p=4,constant=T,irf.periods=20,
                  keep=10000,burnin=5000,VType=1,
                  HP1=0.5,HP2=0.5,HP3=10)
plot(testbvarm,save=F)
IRF(testbvarm,save=F)
forecast(testbvarm,shocks=TRUE,back_data=10,save=FALSE)
}
}

\keyword{BVAR}
\keyword{DSGE}
\concept{BVAR}
\concept{DSGE}