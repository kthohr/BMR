\name{states}

\alias{states}
\alias{states.DSGEVAR}
\alias{states.EDSGE}

\title{Plot State Variables.}

\description{
  Plot the filtered state variables of an estimated DSGE or DSGE-VAR model.
}

\usage{
states(obj,...)

\method{states}{DSGEVAR}(obj,percentiles=c(.05,.50,.95),
       var_names=NULL,use_mean=FALSE,
       save=FALSE,height=13,width=11,...)
\method{states}{EDSGE}(obj,percentiles=c(.05,.50,.95),
       var_names=NULL,use_mean=FALSE,
       save=FALSE,height=13,width=11,...)
}

\arguments{
  For objects of class `EDSGE' or `DSGEVAR',
  \item{obj}{
  An object of class `EDSGE' or `DSGEVAR'.}
	\item{percentiles}{
	Which percentiles of the distribution to use.}
  \item{var_names}{
  Name labels for the states.}
	\item{use_mean}{
	Whether the user would prefer to use the mean of the forecast distribution rather than the middle value in `percentiles'.}
	\item{save}{
	Whether to save the plots.}
	\item{height}{
	If save=FALSE, use this to set the height of the plot.}
	\item{width}{
	If save=FALSE, use this to set the width of the plot.}
  \item{\dots}{
  Additional arguments (not used).}
}

\details{
  This function will work with estimated DSGE and DSGEVAR models.
}

\value{
The function returns a plot of the states with user-selected percentiles, as well as the values used to create the plot; see the vignette for more details on the values returned.
}

\author{Keith O'Hara}

\examples{
\dontrun{
NKMest <- EDSGE(dsgedata,chains=1,cores=1,
                ObserveMat,initialvals,partomats,
                priorform,priorpars,parbounds,par_names,
                optimMethod=c("Nelder-Mead","CG"),
                optimLower=NULL,optimUpper=NULL,
                optimControl=list(maxit=10000),
                DSGEIRFs=TRUE,irf.periods=40,
                scalepar=0.27,keep=50000,burnin=75000,
                tables=TRUE)
#
var_names=c("Output Gap","Output","Inflation","Natural Int",
           "Nominal Int","Labour Supply","Technology","MonetaryPolicy")
states(NKMest,percentiles=c(0.01,0.50,0.99),var_names=var_names)
}
}

\keyword{states}
\concept{states}