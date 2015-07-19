\name{BVARTVP}

\alias{BVARTVP}
\alias{BVARTVP.default}

\title{BVAR with Time-Varying Coefficients.}

\description{
  Estimate a Bayesian VAR with time-varying coefficients.
}

\usage{
BVARTVP(mydata,timelab=NULL,coefprior=NULL,tau=NULL,p=4,
        irf.periods=20,irf.points=NULL,
        keep=10000,burnin=5000,
        XiBeta=1,XiQ=0.01,gammaQ=NULL,
        XiSigma=1,gammaS=NULL)
}

\arguments{
  \item{mydata}{
	A matrix or data frame containing the data series used for estimation; this should be of size T x m.}
	\item{timelab}{
	This is a numeric vector of length T that provides labels for the observations.}
	\item{coefprior}{
	A numeric vector of length m, matrix of size (m p + 1) x m, or a value `NULL', that contains the prior mean-value of each coefficient. Providing a numeric vector of length m will set a zero prior on all coefficients except the own first-lags, which are set according to the elements in `coefprior'. Setting this input to `NULL' will give a random-walk-in-levels prior. 
  
  Note that, when tau is set to `NULL', this input becomes the initial draw for the sampling algorithm, and starting with an explosive draw might be a bad idea.}
	\item{tau}{
	`tau' is the length of the training-sample prior. If this is set a value other than `NULL', it will replace `coefprior' above with the coefficients from a pre-sampling estimation run. Selecting this option also affects the `XiBeta' choice below.}
	\item{p}{
	The number of lags to include of each variable. The default value is 4.}
	\item{irf.periods}{
	An integer value for the horizon of the impulse response calculations, which must be great than zero. The default value is 20.}
	\item{irf.points}{
	A numeric vector of length (0,T]. If the user supplied a `timelab' list above, then this vector should contain points corresponding to that list. The default of `NULL' will mean that all IRFs, for T - tau, will be computed. The IRFs are stored in a 5 dimensional array of size irf.periods x m x m x length(irf.points) x keep. 
  
  If the number of variables, replications, and/or observations is quite large, then calculating all IRFs will take up a lot of memory.  For example, with an IRF horizon of 20, 3 variables, 200 observations, training sample size of 50, and 50000 post-burn-in replications, we have 1,350,000,000 elements to store.}
	\item{keep}{
	The number of Gibbs sampling replications to keep from the sampling run.}
	\item{burnin}{
	The sample burn-in length for the Gibbs sampler.}
	\item{XiBeta}{
	A numeric vector of length 1 or matrix of size (m p + 1) m x (m p + 1) m that contains the prior covariance of each coefficient for beta_0. The structure of Xi_beta corresponds to vec(beta). 
  
  Note that if tau != NULL, `XiBeta' should be a numeric vector of length 1 (scalar) that scales the OLS estimate of the covariance matrix of beta.}
	\item{XiQ}{
	A numeric vector of length 1 or matrix of size (m p + 1) m x (m p + 1) m that contains the location matrix of the inverse-Wishart prior on Q.}
	\item{gammaQ}{
	A numeric vector of length 1 corresponding to the prior degrees of freedom of the Q matrix. The minimum value is (m p + 1) m + 1, and this is the default value, unless tau != NULL, in which case gamma_S = tau.}
	\item{XiSigma}{
	A numeric vector of length 1 or matrix of size m x m that contains the location matrix of the inverse-Wishart prior on Sigma.}
	\item{gammaS}{
	A numeric vector of length 1 corresponding to the prior degrees of freedom of the error covariance matrix. The minimum value is m+1, and this is the default value.}
}

\details{
  For technical details of the model, see the accompanying vignette.
}

\value{
The function returns an object of class BVARTVP, which contains:

  \item{Beta}{
	A matrix of size (m p + 1) m x (T - tau) containing the posterior mean of the coefficient matrix, beta, in vectorised form, for (tau + 1) : T.}
	\item{BDraws}{
	An array of size (m p + 1) x m x keep x (T - tau) which contains the post burn-in draws of beta.}
	\item{Q}{
	A matrix of size (m p + 1) m x (m p + 1) m containing the posterior mean estimates of the covariance matrix Q.}
	\item{QDraws}{
	An array of size (m p + 1) m x (m p + 1) m x keep which contains the post burn-in draws of Q.}
	\item{Sigma}{
	A matrix of size m x m containing the posterior mean estimates of the residual covariance matrix, Sigma.}
	\item{SDraws}{
	An array of size m x m x keep which contains the post burn-in draws of Sigma.}
	\item{IRFs}{
	Let ell =  number of `irf.points' the user selected. `IRFs' is then a five-dimensional object of size irf.periods x m x m x ell x keep containing the impulse response function calculations; the first m refers to responses to the last m shock. 
  
  Note: IRF calculations are based on a one standard deviation shock to each variable.}
	\item{data}{
	The data used for estimation.}
	\item{irf.points}{
	The points in the sample where the user elected to produce IRFs.}
	\item{tau}{
	The length of the training sample.}
}

\author{Keith O'Hara}

\references{
  Koop, Gary and Dimitris Korobilis, ``Bayesian Multivariate Time Series Methods for Empirical Macroeconomics,'' \emph{Mimeo}, 2010.
}

\seealso{
  \code{\link{IRF.BVARTVP}}, \code{\link{plot.BVARTVP}}.
}

\examples{
\dontrun{
irf.points<-c(1979,1988,1997,2006)
yearlab<-seq(1955.00,2010.75,0.25)
USMacroData<-USMacroData[3:226,2:4]
bvartvptest <- BVARTVP(USMacroData,timelab=yearlab,
                coefprior=NULL,tau=80,p=4,
                irf.periods=20,irf.points=irf.points,
                keep=30000,burnin=40000,
                XiBeta=4,XiQ=0.005,gammaQ=NULL,
                XiSigma=1,gammaS=NULL)
}
}

\keyword{BVARTVP}
\concept{BVARTVP}