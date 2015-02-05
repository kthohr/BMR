\name{forecast}

\alias{forecast}
\alias{forecast.BVARM}
\alias{forecast.BVARS}
\alias{forecast.BVARW}
\alias{forecast.CVAR}
\alias{forecast.DSGEVAR}
\alias{forecast.EDSGE}

\title{Forecasting with VAR, DSGE, and DSGE-VAR Models.}

\description{
  Compute forecasts with estimated BVAR, classical VAR, DSGE, and DSGE-VAR models.
}

\usage{
forecast(obj,periods=20,shocks=TRUE,plot=TRUE,
         percentiles=c(.05,.50,.95),useMean=FALSE,
         backdata=0,save=FALSE,height=13,width=11)
forecast(obj,periods=20,plot=TRUE,confint=0.95,
         backdata=0,save=FALSE,height=13,width=11)
forecast(obj,periods=20,plot=TRUE,
         percentiles=c(.05,.50,.95),useMean=FALSE,
         backdata=0,save=FALSE,height=13,width=11)
}

\arguments{
  For objects of class `BVARM', `BVARS', `BVARW', or `DSGEVAR',
  \item{obj}{
  An object of class `BVARM', `BVARS', `BVARW', or `DSGEVAR'.}
	\item{periods}{
	The forecast horizon.}
  \item{shocks}{
	Whether to include uncertainty about future shocks when calculating the distribution of forecasts.}
	\item{plot}{
	Whether to plot the forecasts.}
	\item{percentiles}{
	The percentiles of the distribution the user wants to use.}
	\item{useMean}{
	Whether the user would prefer to use the mean of the forecast distribution rather than the middle value in `percentiles'.}
	\item{backdata}{
	How many 'real' data points to plot before plotting the forecast. A broken line will indicate whether the `real' data ends and the forecast begins.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}\cr
  
  For objects of class `CVAR',
  \item{obj}{
	An object of class `CVAR'.}
	\item{periods}{
	The forecast horizon.}
	\item{plot}{
	Whether to plot the forecasts.}
	\item{confint}{
	The confidence interval to use.}
	\item{backdata}{
	How many 'real' data points to plot before plotting the forecast. A broken line will indicate whether the `real' data ends and the forecast begins.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}\cr
  
  For objects of class `EDSGE',
  \item{obj}{
  An object of class `EDSGE'.}
  \item{periods}{
	The forecast horizon.}
	\item{plot}{
	Whether to plot the forecasts.}
	\item{percentiles}{
	The percentiles of the distribution the user wants to use.}
	\item{useMean}{
	Whether the user would prefer to use the mean of the forecast distribution rather than the middle value in `percentiles'.}
	\item{backdata}{
	How many 'real' data points to plot before plotting the forecast. A broken line will indicate whether the `real' data ends and the forecast begins.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=TRUE, use this to set the height of the plot.}
	\item{width}{
	If save=TRUE, use this to set the width of the plot.}
}

\details{
  This function will work with DSGE, DSGE-VAR, and VAR model estimated with BMR.
}

\value{
The function returns a plot of the forecast with user-selected percentiles, as well as the values used to create the plot; see the vignette for more details on the values returned.
}

\author{Keith O'Hara}

\examples{
\dontrun{
data(BMRVARData)
testbvarw <- BVARW(USMacroData[1:203,2:4],c(0.9,0.95,0.95),p=4,constant=T,
                    irf.periods=20,keep=10000,burnin=5000,
                    XiBeta=4,XiSigma=1,gamma=4)
forecast(testbvarw,periods=10,shocks=T,plot=T,percentiles=c(.05,.50,.95),backdata=10,save=FALSE)
}
}

\keyword{forecast}
\concept{forecast}