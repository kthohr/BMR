\name{modecheck}

\alias{modecheck}

\title{Check the Posterior Mode.}

\description{
  Plot the log posterior around the posterior mode values.
}

\usage{
modecheck(obj,gridsize=1000,scalepar=NULL,plottransform=FALSE,
          save=FALSE,height=13,width=13)
}

\arguments{
  \item{obj}{
	An object of class `EDSGE' or `DSGEVAR'.}
  \item{gridsize}{
  The number of grid points to use when calculating the log posterior around the mode values.}
	\item{scalepar}{
	A value to replace the scale parameter from estimation (`c') when plotting the log posterior.}
	\item{plottransform}{
	Whether to plot the transformed values (i.e., such that the support is unbounded), or to plot the untransformed values.}
	\item{save}{
	Whether to save the plot(s).}
	\item{height}{
	If save=TRUE, use this to set the height of the plot(s).}
	\item{width}{
	If save=TRUE, use this to set the width of the plot(s).}
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
modecheck(NKMest,1000,1,FALSE,save=FALSE)
}
}

\keyword{modecheck}
\concept{modecheck}