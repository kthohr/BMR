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
data(BMRVARData)
testbvarw <- BVARW(USMacroData[,2:4],c(0.9,0.95,0.95),p=4,constant=T,
                    irf.periods=20,keep=10000,burnin=5000,
                    XiBeta=4,XiSigma=1,gamma=4)
IRF(testbvarw,percentiles=c(0.05,0.5,0.95),save=FALSE)
}
}

\keyword{IRF}
\concept{IRF}