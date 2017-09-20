\name{forecast}

\alias{forecast}
\alias{forecast.bvarm}
\alias{forecast.bvars}
\alias{forecast.bvarw}
\alias{forecast.cvar}
\alias{forecast.Rcpp_bvarm}
\alias{forecast.Rcpp_bvars}
\alias{forecast.Rcpp_bvarw}
\alias{forecast.Rcpp_cvar}

\alias{forecast.Rcpp_dsge_gensys}
\alias{forecast.Rcpp_dsge_uhlig}
\alias{forecast.Rcpp_dsgevar_gensys}
\alias{forecast.Rcpp_dsgevar_uhlig}

\title{Forecasting with VAR, DSGE, and DSGE-VAR Models.}

\description{
  Compute forecasts using VAR, DSGE, and DSGE-VAR models.
}

\usage{
forecast(obj,...)

\method{forecast}{Rcpp_bvarm}(obj,periods=20,shocks=TRUE,plot=TRUE,
         var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
\method{forecast}{Rcpp_bvars}(obj,periods=20,shocks=TRUE,plot=TRUE,
         var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
\method{forecast}{Rcpp_bvarw}(obj,periods=20,shocks=TRUE,plot=TRUE,
         var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
\method{forecast}{Rcpp_cvar}(obj,periods=20,shocks=TRUE,plot=TRUE,
         var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)

\method{forecast}{Rcpp_dsge_gensys}(obj,periods=20,shocks=TRUE,plot=TRUE,
         var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
\method{forecast}{Rcpp_dsge_uhlig}(obj,periods=20,shocks=TRUE,plot=TRUE,
         var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)

\method{forecast}{Rcpp_dsgevar_gensys}(obj,periods=20,shocks=TRUE,plot=TRUE,
         var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
\method{forecast}{Rcpp_dsgevar_uhlig}(obj,periods=20,shocks=TRUE,plot=TRUE,
         var_names=NULL,percentiles=c(.05,.50,.95),use_mean=FALSE,
         back_data=0,save=FALSE,height=13,width=11,...)
}

\arguments{
  \item{obj}{
  An object of one of the above classes.}
	\item{periods}{
	The forecast horizon.}
  \item{shocks}{
	Whether to include uncertainty about future shocks when calculating the forecasts.}
	\item{plot}{
	Whether to plot the forecasts.}
    \item{var_names}{
	Variable names.}
	\item{percentiles}{
	The percentiles of the conditional posterior distribution of forecasts to use for plotting.}
	\item{use_mean}{
	Whether to use the mean of the forecast distribution rather than the middle value in `percentiles'.}
	\item{back_data}{
	How many 'real' data points to plot before plotting the forecast. A broken line will indicate whether the `real' data ends and the forecast begins.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=FALSE, use this to set the height of the plot.}
	\item{width}{
	If save=FALSE, use this to set the width of the plot.}\cr
  
  \item{\dots}{
  Additional arguments (not used).}
}

\value{
The function returns a plot of the forecast with user-selected percentiles, as well as the values used to create the plot; see the vignette for more details on the values returned.
}

\author{Keith O'Hara}

\examples{
\dontrun{
data(BMRVARData)
bvar_data <- data.matrix(USMacroData[,2:4])

#

coef_prior <- c(0.9,0.9,0.9)
XiBeta <- 4
XiSigma <- 1
gamma = 4

bvar_obj <- new(bvarw)

#

bvar_obj$build(bvar_data,TRUE,4)
bvar_obj$prior(coef_prior,XiBeta,XiSigma,gamma)
bvar_obj$gibbs(10000,5000)

forecast(bvar_obj,shocks=TRUE,var_names=colnames(bvar_data),back_data=10,save=FALSE)
}
}

\keyword{forecast}
\concept{forecast}