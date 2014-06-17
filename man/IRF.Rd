\name{IRF}

\alias{IRF}
\alias{IRF.BVARM}
\alias{IRF.BVARS}
\alias{IRF.BVARW}
\alias{IRF.CVAR}

\title{Plotting the IRFs.}

\description{
  A function for plotting impulse response funtions using ggplot2.
}

\usage{
IRF(obj,percentiles=c(.05,.50,.95),save=TRUE,height=13,width=13)
}

\arguments{
  \item{obj}{
  An object of class `BVARM', `BVARS', `BVARW', or `CVAR'.}
	\item{percentiles}{
	The percentiles of the distribution the user wants to use.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}
}

\details{
  This function will work with any estimated VAR model from BMR.
}

\value{
The function returns a plot of the IRFs with user-selected percentiles.
}

\author{Keith O'Hara}

\examples{
\dontrun{
IRF(testbvarw,save=TRUE)
}
}

\keyword{IRF}
\concept{IRF}