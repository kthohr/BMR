\name{statespace}

\alias{statespace}
\alias{statespace.gensys}
\alias{statespace.uhlig}
\alias{statespace.SDSGE}

\title{DSGE State-Space Format.}

\description{
  State-space representation of a DSGE model.
}

\usage{
statespace(obj)
}

\arguments{
  \item{obj}{
	An object of type `gensys', `uhlig', or `SDSGE'.}
}

\details{
For technical details, see the accompanying vignette.
}

\value{
The function constructs a state-space representation for the various DSGE solvers in BMR.
}

\author{Keith O'Hara}

\seealso{
  \code{\link{IRF}}, \code{\link{DSGESim}}.
}

\examples{
\dontrun{
mats <- list()
mats$Gamma0 <- Gamma0; mats$Gamma1 <- Gamma1; mats$C <- C; mats$Psi <- Psi; mats$Pi <- Pi 
dsgetest <- SDSGE(mats)
ssmats <- statespace(dsgetest)
}
}

\keyword{statespace}
\concept{statespace}