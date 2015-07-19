\name{gtsplot}

\alias{gtsplot}

\title{Time-series plot using ggplot2.}

\description{
  A simple wrapper function for plotting time-series data with ggplot2.
}

\usage{
gtsplot(X,dates=NULL,rowdates=FALSE,dates.format,
        save=FALSE,height=13,width=11)
}

\arguments{
  \item{X}{
	A matrix or data frame of size T x m containing the relevant time-series data, where m is the number of series.}
	\item{dates}{
	A T x 1 date or character vector containing the relevant date stamps for the data.}
	\item{rowdates}{
	A TRUE or FALSE statement indicating whether the row names of the X matrix contain the date stamps for the data.}
	\item{dates.format}{
	If `dates' is not set to NULL, then indicate what format the dates are in, such as Year-Month-Day.}
  \item{save}{
  Whether to save the plot(s).}
  \item{height}{
  The height of the saved plot(s).}
  \item{width}{
  The width of the saved plot(s).}
}

\value{
  Returns a plot of the data.
}

\author{Keith O'Hara}

\examples{
\dontrun{
data(BMRVARData)
gtsplot(USMacroData[,2:4],dates=USMacroData[,1])
}
}

\keyword{gtsplot}
\concept{gtsplot}