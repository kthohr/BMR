\name{SDSGE}

\alias{SDSGE}
\alias{SDSGE.default}

\title{Solve a DSGE Model.}

\description{
  Solve a DSGE model using one of the available solvers in BMR.
}

\usage{
SDSGE(mats,type=NULL)
}

\arguments{
  \item{mats}{
	The relevant matrices.}
  \item{type}{
  Which solver to use: 1 (gensys) or 2 (uhlig). If type=NULL, the solver will attempt to detect which method to use based on the names of the objects contained in `mats'.}
}

\details{
For technical details of the available solution methods, see the accompanying vignette.
}

\value{
The function returns an object of class `SDSGE', which contains the output of either `gensys' or `uhlig'.
}

\author{Keith O'Hara}

\references{
  Sims, Christopher A., \emph{Solving Linear Rational Expectations Models}, Computational Economics, 20(2), 2002.
  
  Uhlig, Harald, \emph{A Toolkit for Analysing Nonlinear Dynamic Stochastic Models Easily} Computational Methods for the Study of Dynamic Economics, Oxford University Press, 1999.
}

\seealso{
  \code{\link{IRF.SDSGE}}, \code{\link{DSGESim}}.
}

\examples{
\dontrun{
mats <- list()
mats$Gamma0 <- Gamma0; mats$Gamma1 <- Gamma1; mats$C <- C; mats$Psi <- Psi; mats$Pi <- Pi 
dsgetest <- SDSGE(mats)
}
}

\keyword{SDSGE}
\concept{SDSGE}