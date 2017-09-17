\name{uhlig}

\alias{uhlig}
\alias{uhlig.default}

\title{Uhlig's Method.}

\description{
  Solve a DSGE model using Uhlig's method.
}

\usage{
uhlig(A,B,C,D,F,G,H,J,K,L,M,N,whichEig=NULL)
}

\arguments{
  \item{A,B,C,D}{
	The `uhlig' function requires the three blocks of matrices, with 12 matrices in total. The A, B, C, and D matrices form the deterministic block.}
	\item{F,G,H,J,K,L,M}{
	The F, G, and H matrices form the expectational block for the control variables. The J and K matrices are for the `jump' variables, and L and M are for the exogenous shocks.}
	\item{N}{
	The N matrix defines the autoregressive structure of any exogenous shocks.}
	\item{whichEig}{
	The function will return the eigenvalues and (right) eigenvectors used to construct the solution matrices, with the eigenvalues sorted in order of smallest to largest (in absolute value). By default, BMR will select the first (smallest) m eigenvalues (out of a total of 2m eigenvalues). However, if you prefer to select the eigenvalues yourself, then enter a numeric vector of length m indicating which elements of the eigenvalue matrix you wish to use.}
}

\details{
For the technical details of Uhlig's method, see the accompanying vignette.
}

\value{
The function returns an object of class `uhlig', which contains:

  \item{N}{
	The user-specified N matrix, defining the autoregressive structure of any exogenous shocks.}
	\item{P}{
	The P matrix from Uhlig's solution.}
	\item{Q}{
	The Q matrix from Uhlig's solution.}
	\item{R}{
	The R matrix from Uhlig's solution.}
	\item{S}{
	The S matrix from Uhlig's solution.}
	\item{EigenValues}{
	The sorted eigenvalues that form the solution to the P matrix. If a situation of plus/minus infinity in the real part of an eigenvalue (with a corresponding NaN-valued imaginary part) arises, the eigenvalue will be set to 1E+07 +0i.}
	\item{EigenVectors}{
	The eigenvectors corresponding to the sorted eigenvalues.}
}

\author{Keith O'Hara}

\references{
  Uhlig, Harald, \emph{A Toolkit for Analysing Nonlinear Dynamic Stochastic Models Easily} Computational Methods for the Study of Dynamic Economics, Oxford University Press, 1999.
}

\seealso{
  \code{\link{IRF.uhlig}}, \code{\link{DSGESim}}.
}

\examples{
\dontrun{
dsgetest <- uhlig(A,B,C,D,F,G,H,J,K,L,M,N)
}
}

\keyword{uhlig}
\concept{uhlig}