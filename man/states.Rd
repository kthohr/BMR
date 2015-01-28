\name{states}

\alias{states}
\alias{states.DSGEVAR}
\alias{states.EDSGE}

\title{State Plot.}

\description{
  Plot the filtered values of the state variables of an estimated DSGE model.
}

\usage{
states(obj,percentiles=c(.05,.50,.95),
       varnames=NULL,useMean=FALSE,
       save=FALSE,height=13,width=11)
}

\arguments{
  For objects of class `EDSGE',
  \item{obj}{
  An object of class `EDSGE'.}
	\item{percentiles}{
	The percentiles of the distribution the user wants to use.}
  \item{varnames}{
  Name labels for the states.}
	\item{useMean}{
	Whether the user would prefer to use the mean of the forecast distribution rather than the middle value in `percentiles'.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}
}

\details{
  This function will work with DSGE models estimated with BMR.
}

\value{
The function returns a plot of the forecast with user-selected percentiles, as well as the values used to create the plot; see the vignette for more details on the values returned.
}

\author{Keith O'Hara}

\examples{
\dontrun{
states(NKMest)
}
}

\keyword{states}
\concept{states}