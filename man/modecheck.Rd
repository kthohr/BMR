\name{modecheck}

\alias{modecheck}

\title{Check the Posterior Mode.}

\description{
  Plot the log posterior at the posterior mode values.
}

\usage{
modecheck(obj,scalepar=NULL,plottransform=TRUE,save=FALSE,height=13,width=13)
}

\arguments{
  \item{obj}{
	An object of class `EDSGE'.}
	\item{scalepar}{
	A value to replace the scale parameter from estimation (`c') when plotting the log posterior.}
	\item{plottransform}{
	Whether to plot the transformed values (i.e., such that the support is unbounded), or to plot the untransformed values.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}
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
modecheck(NKMest,1,FALSE,save=FALSE)
}
}

\keyword{modecheck}
\concept{modecheck}