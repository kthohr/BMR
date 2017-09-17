\name{bvars}

\alias{bvars}
\alias{Rcpp_bvars}

\title{BVAR with Steady-State Prior.}

\description{
  Estimate a Bayesian VAR with steady-state prior.
}

\usage{
bvar_obj <- new(bvars)
}

\details{
  For technical details of the model, see the accompanying vignette.
}

\author{Keith O'Hara}

\references{
  Villani, Mattias, ``Steady-State Priors for Vector Autoregressions,'' \emph{Journal of Applied Econometrics, 2009, 24} (4), 630--650.
}

\seealso{
  \code{\link{forecast.bvars}}, \code{\link{IRF.bvars}}, \code{\link{plot.bvars}}.
}

\examples{
\dontrun{
data(BMRVARData)
bvar_data <- data.matrix(USMacroData[,2:4])

#

coef_prior <- c(0.9,0.9,0.9)
psi_prior <- matrix(c(3,6,5),nrow=1)

HP1 <- 0.5
HP2 <- 0.5
XiPsi <- 1
gamma = 4

bvar_obj <- new(bvars)

#

bvar_obj$build(bvar_data,TRUE,4)
bvar_obj$prior(coef_prior,HP1,HP2,psi_prior,XiPsi,gamma)
bvar_obj$gibbs(10000,5000)

IRF(bvar_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(bvar_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(bvar_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)
}
}

\keyword{BVAR}
\concept{BVAR}