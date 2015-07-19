\name{prior}

\alias{prior}

\title{Parameterize the Prior Distributions.}

\description{
  Plot and print the moments of a specified prior distribution.
}

\usage{
prior(priorform,priorpars,parname=NULL,moments=TRUE,NR=NULL,NC=NULL)
}

\arguments{
  \item{priorform}{
	This should be a valid prior form for the EDSGE or DSGEVAR functions, such as ``Gamma'' or ``Beta''.}
	\item{priorpars}{
	The relevant parameters of the distribution.}
	\item{parname}{
	A title for the plot.}
	\item{moments}{
  Whether to print the mean, mode, and variance of the distribution.}
  \item{NR}{
  For use with multiple plots. See the vignette for an example.}
  \item{NC}{
  For use with multiple plots. See the vignette for an example.}
}

\details{
  This function can be used when selecting appropriate prior distributions for the EDSGE and DSGEVAR functions.
}

\author{Keith O'Hara}

\examples{
\dontrun{
prior("Normal",c(0,1))
prior("Gamma",c(2,2))
}
}

\keyword{prior}
\concept{prior}