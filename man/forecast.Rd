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
  Compute forecasts using VAR, DSGE, and DSGE-VAR models.
}

\usage{
forecast(obj,...)

\method{forecast}{BVARM}(obj,periods=20,shocks=TRUE,plot=TRUE,
         percentiles=c(.05,.50,.95),useMean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
\method{forecast}{BVARS}(obj,periods=20,shocks=TRUE,plot=TRUE,
         percentiles=c(.05,.50,.95),useMean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
\method{forecast}{BVARW}(obj,periods=20,shocks=TRUE,plot=TRUE,
         percentiles=c(.05,.50,.95),useMean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
\method{forecast}{DSGEVAR}(obj,periods=20,shocks=TRUE,plot=TRUE,
         percentiles=c(.05,.50,.95),useMean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)

\method{forecast}{CVAR}(obj,periods=20,plot=TRUE,confint=0.95,
         back_data=0,save=FALSE,height=13,width=11,...)

\method{forecast}{EDSGE}(obj,periods=20,plot=TRUE,
         percentiles=c(.05,.50,.95),useMean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
}

\arguments{
  For objects of class `BVARM', `BVARS', `BVARW', or `DSGEVAR',
  \item{obj}{
  An object of the above class.}
	\item{periods}{
	The forecast horizon.}
  \item{shocks}{
	Whether to include uncertainty about future shocks when calculating the distribution of forecasts.}
	\item{plot}{
	Whether to plot the forecasts.}
	\item{percentiles}{
	The percentiles of the conditional posterior distribution of forecasts to use for plotting.}
	\item{useMean}{
	Whether the user would prefer to use the mean of the forecast distribution rather than the middle value in `percentiles'.}
	\item{back_data}{
	How many 'real' data points to plot before plotting the forecast. A broken line will indicate whether the `real' data ends and the forecast begins.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=FALSE, use this to set the height of the plot.}
	\item{width}{
	If save=FALSE, use this to set the width of the plot.}\cr
  
  For objects of class `CVAR',
	\item{confint}{
	The confidence interval to use.}\cr
  
  \item{\dots}{
  Additional arguments (not used).}
}

\details{
  This function will work with DSGE, DSGE-VAR, and VAR models estimated with BMR.
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
forecast(testbvarw,periods=10,shocks=T,plot=T,percentiles=c(.05,.50,.95),back_data=10,save=FALSE)
}
}

\keyword{forecast}
\concept{forecast}