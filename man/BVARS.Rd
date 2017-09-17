\name{bvars}

\alias{bvars}
\alias{bvars.default}

\title{BVAR with Steady-State Prior.}

\description{
  Estimate a Bayesian VAR with steady-state prior.
}

\usage{
bvars(mydata,psiprior=NULL,coefprior=NULL,p=4,
      irf.periods=20,keep=10000,burnin=1000,
      XiPsi=1,HP1=0.5,HP4=2,gamma=NULL)
}

\arguments{
  \item{mydata}{
  A matrix or data frame containing the data series used for estimation; this should be of size T x m.}
  \item{psiprior}{
  A numeric vector of length m that contains the prior mean of each series found in `mydata'. The user MUST specify this prior, or the function will return an error.}
	\item{coefprior}{
	A numeric vector of length m, matrix of size (m p) x m, or a value of `NULL', that contains the prior mean-value of each coefficient. Providing a numeric vector of length m will set a zero prior on all coefficients except the own first-lags, which are set according to the elements in `coefprior'. Setting this input to `NULL' will give a random-walk-in-levels prior.}
	\item{p}{
	The number of lags to include of each variable. The default value is 4.}
	\item{irf.periods}{
	An integer value for the horizon of the impulse response calculations, which must be greater than zero. The default value is 20.}
	\item{keep}{
	The number of Gibbs sampling replications to keep from the sampling run.}
	\item{burnin}{
	The sample burn-in length for the Gibbs sampler.}
  \item{XiPsi}{
	A numeric vector of length 1 or matrix of size m x m that defines the prior location matrix of Psi.}
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
The function returns an object of class bvars, which contains:

  \item{Beta}{
	A matrix of size (m p) x m containing the posterior mean of the coefficient matrix (beta).}
	\item{BDraws}{
	An array of size (m p) x m x keep which contains the post burn-in draws of beta.}
  \item{Psi}{
  A matrix of size 1 x m containing the posterior mean estimate of the unconditional mean matrix (Psi).}
	\item{PDraws}{
	An array of size 1 x m x keep which contains the post burn-in draws of Psi.}
	\item{Sigma}{
	A matrix of size m x m containing the posterior mean estimate of the residual covariance matrix (Sigma).}
	\item{SDraws}{
	An array of size m x m x keep which contains the post burn-in draws of Sigma.}
	\item{IRFs}{
	A four-dimensional object of size irf.periods x m x keep x m containing the impulse response function calculations; the first m refers to responses to the last m shock. 
  
  Note: IRF calculations are based on a one standard deviation shock to each variable.}
	\item{data}{
	The data used for estimation.}
}

\author{Keith O'Hara}

\references{
  Villani, Mattias, ``Steady-State Priors for Vector Autoregressions,'' \emph{Journal of Applied Econometrics, 2009, 24} (4), 630--650.
}

\seealso{
  \code{\link{forecast.bvars}}, \code{\link{IRF.bvars}}, \code{\link{plot.bvars}}.
}

\examples{
\dontrun{
data(BMRVARData)
mycfp = c(0.9,0.95,0.95)
mypsi = c(3,6,5)
testbvars <- bvars(USMacroData[,2:4],mypsi,mycfp,p=4,irf.periods=20,
                   keep=20000,burnin=5000,XiPsi=1,HP1=0.5,HP4=2,gamma=NULL)
plot(testbvars,save=F)
IRF(testbvars,save=F)
forecast(testbvars,shocks=T,back_data=10,save=F)
}
}

\keyword{BVAR}
\concept{BVAR}