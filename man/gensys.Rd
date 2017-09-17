\name{gensys}

\alias{gensys}
\alias{gensys.default}

\title{Gensys Solver.}

\description{
  Solve a DSGE model using Sims' method.
}

\usage{
gensys(Gamma0,Gamma1,C,Psi,Pi)
}

\arguments{
  \item{Gamma0}{
	Coefficients on present-time variables.}
	\item{Gamma1}{
	Coefficients on lagged variables.}
	\item{C}{
	Intercept terms.}
	\item{Psi}{
  Coefficients on any exogenous shocks.}
  \item{Pi}{
  One-step-ahead expectational errors.}
}

\details{
For the technical details of Sims' method, see the accompanying vignette.
}

\value{
The function returns an object of class `gensys', which contains:

  \item{G1}{
	Autoregressive solution matrix.}
	\item{Cons}{
	Intercept terms.}
	\item{impact}{
	Coefficients on the exogenous shocks.}
	\item{eu}{
	A 2 x 1 vector indicating existence and uniqueness (respectively) of the solution. A value of 1 can be read as `yes', while 0 is `no'.}
	\item{Psi}{
	User-specified shock matrix.}
	\item{Pi}{
	User-specified expectational errors matrix.}
}

\author{Keith O'Hara}

\references{
  Sims, Christopher A., \emph{Solving Linear Rational Expectations Models}, Computational Economics, 20(2), 2002.
}

\seealso{
  \code{\link{IRF.gensys}}, \code{\link{DSGESim}}.
}

\examples{
\dontrun{
dsgetest <- gensys(Gamma0,Gamma1,C,Psi,Pi)
}
}

\keyword{gensys}
\concept{gensys}