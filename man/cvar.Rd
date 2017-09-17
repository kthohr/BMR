\name{cvar}

\alias{cvar}
\alias{cvar.default}

\title{Classical VAR.}

\description{
  OLS estimation of a VAR model with bootstrapped IRFs.
}

\usage{
var_obj <- new(cvar)
}

\details{
  For technical details of the model, see the accompanying vignette.
}

\author{Keith O'Hara}

\seealso{
  \code{\link{forecast.cvar}}, \code{\link{IRF.cvar}}.
}

\examples{
\dontrun{
data(BMRVARData)
var_data <- data.matrix(USMacroData[,2:4])

var_obj = new(cvar)

#

var_obj$build(var_data,TRUE,4)
var_obj$estim()
var_obj$boot(10000)

IRF(var_obj,20,var_names=colnames(USMacroData),save=FALSE)
plot(var_obj,var_names=colnames(USMacroData),save=FALSE)
forecast(var_obj,shocks=TRUE,var_names=colnames(USMacroData),back_data=10,save=FALSE)
}
}

\keyword{VAR}
\concept{VAR}