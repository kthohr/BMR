\name{gpacf}

\alias{gpacf}

\title{Partial Autocorrelation Function.}

\description{
  Plotting the partial autocorrelation function using ggplot2.
}

\usage{
gpacf(y,lags=10,ci=.95,plot=T,barcolor="darkred",
      names=FALSE,save=FALSE,height=12,width=12)
}

\arguments{
  \item{y}{
  A matrix or data frame of size T x m containing the relevant series.}
	\item{lags}{
	The number of lags to plot.}
	\item{ci}{
	A numeric value betweein 0 and 1 specifying the confidence interval to use; the default value is 0.95.}
	\item{barcolor}{
	The color of the bars.}
  \item{names}{
	Whether to plot the names of the series.}
	\item{save}{
	Whether to save the plots. The default is `FALSE'.}
	\item{height}{
	If save = TRUE, use this to set the height of the plot.}
	\item{width}{
	If save = TRUE, use this to set the width of the plot.}
}

\details{
Plotting the partial autocorrelation function.
}

\value{
  Returns a plot of the PACF.
}

\author{Keith O'Hara}

\examples{
\dontrun{
data(BMRVARData)
gpacf(USMacroData[,2:4],lags=12,ci=0.95,plot=T,barcolor="darkred",
      names=T,save=T,height=6,width=12)
}
}

\keyword{gpacf}
\keyword{PACF}
\concept{ACF}