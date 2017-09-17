\name{modecheck}

\alias{modecheck}
\alias{modecheck.default}
\alias{modecheck.EDSGE}
\alias{modecheck.DSGEVAR}

\title{Check the Posterior Mode.}

\description{
  Plot the log posterior around the posterior mode values.
}

\usage{
modecheck(obj,...)

\method{modecheck}{DSGEVAR}(obj,gridsize=200,scalepar=NULL,plottransform=FALSE,par_names=NULL,
  save=FALSE,height=13,width=13,...)
\method{modecheck}{EDSGE}(obj,gridsize=200,scalepar=NULL,plottransform=FALSE,par_names=NULL,
  save=FALSE,height=13,width=13,...)
}

\arguments{
  \item{obj}{
	An object of class `EDSGE' or `DSGEVAR'.}
  \item{gridsize}{
  The number of grid points to use when calculating the log posterior around the mode values.}
	\item{scalepar}{
	A value to replace the scaling parameter from estimation (`c') when plotting the log posterior.}
	\item{plottransform}{
	Whether to plot the transformed values (i.e., such that the support of each parameter is unbounded), or to plot the untransformed values.}
  \item{par_names}{
  A vector of expressions containing the name of each parameter.}
	\item{save}{
	Whether to save the plot(s).}
	\item{height}{
	If save=FALSE, use this to set the height of the plot(s).}
	\item{width}{
	If save=FALSE, use this to set the width of the plot(s).}
  \item{\dots}{
  Additional arguments (not used).}
}

\details{
  For an example, see the accompanying vignette.
}

\value{
The function will plot the value of the log posterior around the posterior mode.
}

\author{Keith O'Hara}

\examples{
\dontrun{
modecheck(NKMest,200,1,FALSE,save=FALSE)
}
}

\keyword{modecheck}
\concept{modecheck}