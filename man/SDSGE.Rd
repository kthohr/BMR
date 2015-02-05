\name{SDSGE}

\alias{SDSGE}

\title{Solve DSGE.}

\description{
  Solve a DSGE Model using Uhlig's method.
}

\usage{
SDSGE(A,B,C,D,F,G,H,J,K,L,M,N,whichEig=NULL)
}

\arguments{
  \item{A,B,C,D}{
	The SDSGE function requires the three blocks from Uhlig's method, 12 matrices in total. The A, B, C, and D matrices form the deterministic block.}
	\item{F,G,H,J,K,L,M}{
	The F, G, and H matrices form the expectational block for the control variables. The J and K matrices are for the `jump' variables, and L and M are for the exogenous shocks.}
	\item{N}{
	The N matrix defines the autoregressive nature of the exogenous shocks.}
	\item{whichEig}{
	The SDSGE function will return the eigenvalues and (right) eigenvectors used to construct the P matrix, in order of smallest eigenvalues to the largest. By default, BMR will select the first (smallest) m eigenvalues (out of a total of 2m eigenvalues). If you prefer to select the eigenvalues yourself, then enter a numeric vector of length m indicating which elements of the eigenvalue numeric vector you want to use, and BMR will switch to the corresponding eigenvectors automatically.}
}

\details{
For the technical details of Uhlig's method, see the accompanying vignette.
}

\value{
The function returns an object of class SDSGE, which contains:

  \item{N}{
	The user-specified N matrix, defining the autoregressive nature of any exogenous shocks.}
	\item{P}{
	The P matrix from Uhlig's solution.}
	\item{Q}{
	The Q matrix from Uhlig's solution.}
	\item{R}{
	The R matrix from Uhlig's solution.}
	\item{S}{
	The S matrix from Uhlig's solution.}
	\item{EigenValues}{
	The sorted eigenvalues that form the solution to the P matrix. If a situation of plus/minus infinity in the real part of an eigenvalue (with a corresponding NaN-valued imaginary part) arises, the eigenvalue will be set to 1E+07 +0i. (Armadillo doesn't like sorting vectors with infinite values.)}
	\item{EigenVectors}{
	The eigenvectors that correspond to the sorted eigenvalues.}
}

\author{Keith O'Hara}

\references{
  Uhlig, Harald, \emph{A Toolkit for Analysing Nonlinear Dynamic Stochastic Models Easily} Computational Methods for the Study of Dynamic Economics, Oxford University Press, 1999.
}

\seealso{
  \code{\link{IRF}}, \code{\link{DSGESim}}.
}

\examples{
\dontrun{
dsgetest <- SDSGE(A,B,C,D,F,G,H,J,K,L,M,N)
}
}

\keyword{SDSGE}
\concept{SDSGE}