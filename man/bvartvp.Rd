\name{bvartvp}

\alias{bvartvp}
\alias{Rcpp_bvartvp}

\title{BVAR with Time-Varying Coefficients.}

\description{
  Estimate a Bayesian VAR with time-varying coefficients.
}

\usage{
bvar_obj <- new(bvartvp)
}

\details{
  For technical details of the model, see the accompanying vignette.
}

\author{Keith O'Hara}

\references{
  Koop, Gary and Dimitris Korobilis, ``Bayesian Multivariate Time Series Methods for Empirical Macroeconomics,'' \emph{Mimeo}, 2010.
}

\seealso{
  \code{\link{IRF.bvartvp}}, \code{\link{plot.bvartvp}}.
}

\examples{
\dontrun{
data(BMRVARData)
bvar_data <- data.matrix(USMacroData[,2:4])

#

tau <- 80

XiBeta <- 4
XiQ <- 0.005
gammaQ <- tau
XiSigma <- 1
gammaS = 4

which_irfs = c(17,53,89,125)

bvar_obj <- new(bvartvp)

#

bvar_obj$build(bvar_data,TRUE,1)
bvar_obj$prior(tau,XiBeta,XiQ,gammaQ,XiSigma,gammaS)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,which_irfs,var_names=colnames(USMacroData),save=FALSE)
plot(bvar_obj,var_names=colnames(USMacroData),save=FALSE)

}
}

\keyword{bvartvp}
\concept{bvartvp}