\name{bvarm}

\alias{bvarm}
\alias{Rcpp_bvarm}

\title{BVAR with Minnesota Prior.}

\description{
    A Bayesian VAR with Minnesota prior.
}

\usage{
bvar_obj <- new(bvarm)
}

\details{
  For technical details of the model, see the accompanying vignette.
}

\author{Keith O'Hara}

\references{
  Canova, Fabio, \emph{Methods for Applied Macroeconomic Research,} Princeton, New Jersey: Princeton University Press, 2007.
  
  Koop, Gary and Dimitris Korobilis, ``Bayesian Multivariate Time Series Methods for Empirical Macroeconomics,'' \emph{Mimeo}, 2010.
}

\seealso{
  \code{\link{forecast.bvarm}}, \code{\link{IRF.bvarm}}, \code{\link{plot.bvarm}}.
}

\examples{
\dontrun{
data(BMRVARData)
bvar_data <- data.matrix(USMacroData[,2:4])

#

coef_prior <- c(0.9,0.9,0.9)

bvar_obj <- new(bvarm)

bvar_obj$build(bvar_data,TRUE,4)
bvar_obj$prior(coef_prior,1,1,0.5,0.5,100.0,1.0)
bvar_obj$gibbs(10000)

IRF(bvar_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(bvar_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)
}
}

\keyword{BVAR}
\concept{BVAR}
