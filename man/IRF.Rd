\name{IRF}

\alias{IRF}

\alias{IRF.bvarm}
\alias{IRF.bvars}
\alias{IRF.bvartvp}
\alias{IRF.bvarw}
\alias{IRF.cvar}

\alias{IRF.DSGEVAR}
\alias{IRF.EDSGE}

\alias{IRF.gensys}
\alias{IRF.uhlig}
\alias{IRF.SDSGE}

\title{Plotting IRFs.}

\description{
  A function for plotting impulse response funtions using ggplot2.
}

\usage{
  IRF(obj,...)

  \method{IRF}{bvarm}(obj,percentiles=c(.05,.50,.95),which_shock=NULL,which_response=NULL,save=FALSE,height=13,width=13,...)
  \method{IRF}{bvars}(obj,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13,...)
  \method{IRF}{bvarw}(obj,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13,...)
  \method{IRF}{cvar}(obj,percentiles=c(.05,.50,.95),save=FALSE,height=13,width=13,...)

  \method{IRF}{bvartvp}(obj,whichirfs=NULL,percentiles=c(.05,.50,.95),
  save=FALSE,height=13,width=13,...)
    
  \method{IRF}{DSGEVAR}(obj,var_names=NULL,percentiles=c(.05,.50,.95),comparison=TRUE,
    save=FALSE,height=13,width=13,...)
  \method{IRF}{EDSGE}(obj,observableIRFs=FALSE,var_names=NULL,percentiles=c(.05,.50,.95),
    save=FALSE,height=13,width=13,...)

  \method{IRF}{gensys}(obj,shocks,irf.periods=20,var_names=NULL,
    plot=TRUE,save=FALSE,height=13,width=13,...)
  \method{IRF}{uhlig}(obj,shocks,irf.periods=20,var_names=NULL,
    plot=TRUE,save=FALSE,height=13,width=13,...)
  \method{IRF}{SDSGE}(obj,shocks,irf.periods=20,var_names=NULL,
    plot=TRUE,save=FALSE,height=13,width=13,...)
}

\arguments{
  For objects of class `bvarm', `bvars', `bvarw', or `cvar',
  \item{obj}{
  An object of the above class.}
	\item{percentiles}{
	The percentiles of the distribution the user wants to use.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=FALSE, use this to set the height of the plot.}
	\item{width}{
	If save=FALSE, use this to set the width of the plot.}\cr
  
  For objects of class `bvartvp',
  \item{whichirfs}{
  Which IRFs to plot. (The default is to plot all of the IRFs contained in the estimation object.)}\cr
  
  For objects of class `DSGEVAR',
  \item{var_names}{
  A character vector with the names of the relevant variables.}
  \item{comparison}{
  Whether to plot corresponding DSGE model IRFs.}\cr
  
  For objects of class `EDSGE',
  \item{observableIRFs}{
  Whether to plot the IRFs relating to the state variables, or the implied IRFs of the observable series.}\cr
  
  For objects of class `gensys', `uhlig', or `SDSGE',
  \item{shocks}{
	A numeric vector containing the standard deviations of the shocks.}
  \item{irf.periods}{
  The horizon of the IRFs.}
  \item{plot}{
  Whether to plot the IRFs.}\cr
  
  \item{\dots}{
  Additional arguments (not used).}
}

\details{
  This function will work with any estimated VAR model, estimated DSGE or DSGE-VAR model, or solved DSGE model.
}

\value{
The function returns a plot of the IRFs with user-selected percentiles.
}

\author{Keith O'Hara}

\examples{
\dontrun{
data(BMRVARData)
testbvarw <- bvarw(USMacroData[,2:4],c(0.9,0.95,0.95),p=4,constant=T,
                    irf.periods=20,keep=10000,burnin=5000,
                    XiBeta=4,XiSigma=1,gamma=4)
IRF(testbvarw,percentiles=c(0.05,0.5,0.95),save=FALSE)
}
}

\keyword{IRF}
\concept{IRF}