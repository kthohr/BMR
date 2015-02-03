\name{IRF}

\alias{IRF}
\alias{IRF.BVARM}
\alias{IRF.BVARS}
\alias{IRF.BVARTVP}
\alias{IRF.BVARW}
\alias{IRF.CVAR}
\alias{IRF.DSGEVAR}
\alias{IRF.EDSGE}
\alias{IRF.SDSGE}


\title{Plotting IRFs.}

\description{
  A function for plotting impulse response funtions using ggplot2.
}

\usage{
IRF(obj,percentiles=c(.05,.50,.95),save=TRUE,height=13,width=13)
IRF(obj,whichirfs=NULL,percentiles=c(.05,.50,.95),
    save=FALSE,height=13,width=13)
IRF(obj,varnames=NULL,percentiles=c(.05,.50,.95),comparison=TRUE,
    save=TRUE,height=13,width=13)
IRF(obj,ObservableIRFs=TRUE,varnames=NULL,percentiles=c(.05,.50,.95),
    save=TRUE,height=13,width=13)
IRF(obj,shocks,irf.periods=20,varnames=NULL,plot=TRUE,
    save=FALSE,height=13,width=13)
}

\arguments{
  For objects of class `BVARM', `BVARS', `BVARW', or `CVAR',
  \item{obj}{
  An object of class `BVARM', `BVARS', `BVARW', or `CVAR'.}
	\item{percentiles}{
	The percentiles of the distribution the user wants to use.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}\cr
  
  For objects of class `BVARTVP',
  \item{obj}{
  An object of class `BVARTVP'.}
  \item{whichirfs}{
  Which IRFs to plot. The default is all of those which the user chose to calculate during estimation.}
  \item{percentiles}{
	The percentiles of the distribution the user wants to use.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}\cr
  
  For objects of class `DSGEVAR',
  \item{obj}{
  An object of class `DSGEVAR'.}
  \item{varnames}{
  A character vector with the names of the relevant variables.}
  \item{percentiles}{
	The percentiles of the distribution the user wants to use.}
  \item{comparison}{
  Whether to plot corresponding DSGE model IRFs.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}\cr
  
  For objects of class `EDSGE',
  \item{obj}{
  An object of class `EDGSE'.}
  \item{ObservableIRFs}{
  Whether to plot the IRFs relating to the state variables, or the implied IRFs of the observable series.}
  \item{varnames}{
  A character vector with the names of the relevant variables.}
  \item{percentiles}{
	The percentiles of the distribution the user wants to use.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}\cr
  
  For objects of class `SDSGE',
  \item{obj}{
  An object of class `SDSGE'.}
  \item{shocks}{
	A numeric vector containing the standard deviations of the shocks.}
  \item{irf.periods}{
  The horizon of the IRFs.}
  \item{varnames}{
  A character vector with the names of the relevant variables.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}\cr
}

\details{
  This function will work with any estimated VAR model, estimated DSGE or DSGEVAR model, or solved DSGE model.
}

\value{
The function returns a plot of the IRFs with user-selected percentiles.
}

\author{Keith O'Hara}

\examples{
\dontrun{
data(BMRVARData)
testbvarw <- BVARW(USMacroData[,2:4],c(0.9,0.95,0.95),p=4,constant=T,
                    irf.periods=20,keep=10000,burnin=5000,
                    XiBeta=4,XiSigma=1,gamma=4)
IRF(testbvarw,percentiles=c(0.05,0.5,0.95),save=FALSE)
}
}

\keyword{IRF}
\concept{IRF}