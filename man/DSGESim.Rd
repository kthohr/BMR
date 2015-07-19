\name{DSGESim}

\alias{DSGESim}
\alias{DSGESim.gensys}
\alias{DSGESim.uhlig}
\alias{DSGESim.SDSGE}

\title{Simulate from a DSGE Model.}

\description{
  Simulate a data series using a solved DSGE model.
}

\usage{
DSGESim(obj,...)

\method{DSGESim}{gensys}(obj,shocks.cov,sim.periods,burnin=NULL,seedval=1122,hpfiltered=FALSE,lambda=1600,...)
\method{DSGESim}{uhlig}(obj,shocks.cov,sim.periods,burnin=NULL,seedval=1122,hpfiltered=FALSE,lambda=1600,...)
\method{DSGESim}{SDSGE}(obj,shocks.cov,sim.periods,burnin=NULL,seedval=1122,hpfiltered=FALSE,lambda=1600,...)
}

\arguments{
  \item{obj}{
	An object of class `SDSGE', `gensys', or `uhlig'. The user should first solve a model using one of the solver functions (`SDSGE', `gensys', or `uhlig'), then pass the solution to `DSGESim'.}
	\item{shocks.cov}{
	A matrix of size k x k that describes the covariance structure of the model shocks.}
	\item{sim.periods}{
	The number of simulation periods the function should return.}
	\item{burnin}{
	The length of sample burn-in. The default, `burnin = NULL', will set burn-in to one-half of the number given in `sim.periods'.}
  \item{seedval}{
  Seed the random number generator.}
	\item{hpfiltered}{
	Whether to pass the simulated series through a Hodrick-Prescott filter before retuning.}
	\item{lambda}{
	If `hpfiltered = TRUE', this is the value of the smoothing parameter in the H-P filter.}
  \item{\dots}{
  Additional arguments (not used).}
}

\details{
  For an example, see the accompanying vignette.
}

\value{
The function will return a matrix of simulated observations from a solved DSGE model.
}

\author{Keith O'Hara}

\seealso{
  \code{\link{gensys},\link{uhlig},\link{SDSGE}}
}

\examples{
\dontrun{
dsgetest <- gensys(Gamma0,Gamma1,C,Psi,Pi)
dsgetestsim <- DSGESim(dsgetest,Sigma,200,200)
}
}

\keyword{DSGESim}
\concept{DSGESim}