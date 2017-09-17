\name{bvarw}

\alias{bvarw}
\alias{bvarw.default}

\title{BVAR with normal-inverse-Wishart Prior.}

\description{
  Estimate a Bayesian VAR with normal-inverse-Wishart prior.
}

\usage{
bvar_obj <- new(bvarw)
}

\details{
  For technical details of the model, see the accompanying vignette.
}

\author{Keith O'Hara}

\references{
  Koop, Gary and Dimitris Korobilis, ``Bayesian Multivariate Time Series Methods for Empirical Macroeconomics,'' \emph{Mimeo}, 2010.
}

\seealso{
  \code{\link{forecast.bvarw}}, \code{\link{IRF.bvarw}}, \code{\link{plot.bvarw}}.
}

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

IRF(bvar_obj,20,var_names=colnames(bvar_data),save=FALSE)
plot(bvar_obj,var_names=colnames(bvar_data),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(bvar_data),back_data=10,save=FALSE)
}
}

\keyword{BVAR}
\concept{BVAR}