\name{DSGESim}

\alias{DSGESim}
\alias{DSGESim.SDSGE}

\title{Simulate DSGE.}

\description{
  Simulate a solved DSGE model.
}

\usage{
DSGESim(obj,seedval=1122,shocks,sim.periods,burnin=NULL,hpfiltered=FALSE,lambda=1600)
}

\arguments{
  \item{obj}{
	An object of class `SDSGE'. The user should first solve a model using the `SDSGE' function, then pass the solution to the `DSGESim' function.}
  \item{seedval}{
  Seed the random number generator.}
	\item{shocks}{
	A numeric vector of length k or a matrix of size k x k with the standard deviation of each shock.}
	\item{sim.periods}{
	The number of periods the user wishes to simulate for and keep.}
	\item{burnin}{
	The length of sample burn-in. The default, `burnin = NULL', will set this to one-half of the figure given in `sim.periods'.}
	\item{hpfiltered}{
	Whether to pass the simulated series through a Hodrick-Prescott filter before retuning it.}
	\item{lambda}{
	If `hpfiltered = TRUE', this is the value of the smoothing parameter in the H-P filter.}
}

\details{
  For an example, see the accompanying vignette.
}

\value{
The function will return a matrix of simulated observations from a solved DSGE model.
}

\author{Keith O'Hara}

\seealso{
  \code{\link{SDSGE}}
}

\examples{
\dontrun{
dsgetestsim <- DSGESim(dsgetest,1122,1,200,200,hpfiltered=FALSE)
}
}

\keyword{DSGESim}
\concept{DSGESim}