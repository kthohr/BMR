\name{BVARS}

\alias{BVARS}
\alias{plot.BVARS}

\title{BVAR with Steady-State Prior.}

\description{
  Estimating a Bayesian VAR with a steady-state prior.
}

\usage{
BVARS(mydata,psiprior=NULL,coefprior=NULL,p=4,
      irf.periods=20,keep=10000,burnin=1000,
      XiPsi=1,HP1=0.5,HP4=2,gamma=NULL)
\method{plot}{BVARS}(obj,type=1,plotSigma=TRUE,save=TRUE,height=13,width=13)
}

\arguments{
  \item{mydata}{
  A matrix or data frame containing the series to be used in estimation, and should be of size T x m.}
  \item{psiprior}{
  A numeric vector of length m that contains the prior mean of each series found in `mydata'. The user MUST specify this prior, or the function will return an error requesting as such.}
	\item{coefprior}{
	A numeric vector of length m, matrix of size (m p) x m, or logical value `NULL', that contains the prior mean of each coefficient. Only providing a numeric vector of length m will make BMR set a zero prior on all coefficients, except the own first-lags, which are set according to the elements in `coefprior'. Setting this to `NULL' will mean a random walk in levels prior.}
	\item{p}{
	The number of lags to include of each variable. The default value is 4.}
	\item{irf.periods}{
	An integer value for the horizon of the impulse response calculations; this value must be greater than zero. The default value is 20.}
	\item{keep}{
	The number of Gibbs sampling replications to keep from the sampling run.}
	\item{burnin}{
	The sample burn-in period for the Gibbs sampler.}
  \item{XiPsi}{
	A numeric vector of length 1 or matrix of size m x m that defines the prior location matrix of psi.}
	\item{HP1}{
	H_1 from section 3 of the vignette.}
	\item{HP4}{
	H_4 from section 3 of the vignette.}
	\item{gamma}{
	A numeric vector of length 1 corresponding to the prior degrees of freedom of the error covariance matrix. The minimum value is m+1, and this is the default value.}
}

\details{
  For technical details of the model, see the accompanying vignette.
}

\value{
The function returns an object of class BVARS, which contains:

  \item{Beta}{
	A matrix of size (m p) x m containing the mean posterior estimate of the coefficient matrix, beta.}
	\item{BDraws}{
	An array of size (m p) x m x keep which contains the post-burn-in draws of beta.}
  \item{Psi}{
  A matrix of size 1 x m containing the mean posterior estimate of the unconditional mean matrix, Psi.}
	\item{PDraws}{
	An array of size 1 x m x keep which contains the post-burn-in draws of Psi.}
	\item{Sigma}{
	A matrix of size m x m containing the mean posterior estimate of the residual covariance matrix, Sigma.}
	\item{SDraws}{
	An array of size m x m x keep which contains the post-burn-in draws of Sigma.}
	\item{IRFs}{
	A four-dimensional object of size irf.periods x m x keep x m containing the impulse response function calculations; the first m refers to the responses to the last m shock. Note: IRF calculations are based on a one standard deviation shock to each variable.}
	\item{data}{
	The data used for estimation.}
}

\author{Keith O'Hara}

\references{
  Villani, Mattias, ``Steady-State Priors for Vector Autoregressions,'' \emph{Journal of Applied Econometrics, 2009, 24} (4), 630–650.
}

\seealso{
  \code{\link{plot}}
}

\examples{
\dontrun{
data(BMRVARData)
mycfp = c(0.9,0.95,0.95)
mypsi = c(3,6,5)
testbvars <- BVARS(USMacroData[,2:4],mypsi,mycfp,p=4,irf.periods=20,keep=20000,burnin=5000,XiPsi=1,HP1=0.5,HP4=2,gamma=NULL)
}
}

\keyword{BVAR}
\concept{BVAR}