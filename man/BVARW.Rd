\name{BVARW}

\alias{BVARW}
\alias{plot.BVARW}

\title{BVAR with normal-inverse-Wishart Prior.}

\description{
  Estimate a Bayesian VAR with normal-inverse-Wishart prior.
}

\usage{
BVARW(mydata,cores=1,coefprior=NULL,p=4,constant=TRUE,
      irf.periods=20,keep=10000,burnin=1000,
      XiBeta=1,XiSigma=1,gamma=NULL)
\method{plot}{BVARW}(obj,type=1,plotSigma=TRUE,save=TRUE,height=13,width=13)
}

\arguments{
  \item{mydata}{
  A matrix or data frame containing the series to be used in estimation, and should be of size T x m.}
  \item{cores}{
  A positive integer value indicating the number of CPU cores that should be used in estimation. DO NOT ENTER MORE CORES THAN YOUR COMPUTER CAN HANDLE! If in doubt, set cores = 1, which is the default.}
	\item{coefprior}{
	A numeric vector of length m, matrix of size (m p + 1_c) x m, or a logical value of `NULL', that contains the prior mean of each coefficient. Only providing a numeric vector of length m will make BMR set a zero prior on all coefficients except the own first-lags, which are set according to the elements in `coefprior'. Setting this to `NULL' will mean a random walk in levels prior.}
	\item{p}{
	The number of lags to include of each variable. The default value is 4.}
	\item{constant}{
	A logical statement on whether to include a constant vector in the model. The default is `TRUE', and the alternative is `FALSE'.}
	\item{irf.periods}{
	An integer value for the horizon of the impulse response calculations; this value must be greater than zero. The default value is 20.}
	\item{keep}{
	The number of Gibbs sampling replications to keep from the sampling run.}
	\item{burnin}{
	The sample burn-in period for the Gibbs sampler.}
	\item{XiBeta}{
	A numeric vector of length 1 or matrix of size (m p + 1_c) m x (m p + 1_c) m comprising the prior covariance of each coefficient. The structure of Xi_beta corresponds to vec(beta).}
	\item{XiSigma}{
	A numeric vector of length 1 or matrix of size m x m that contains the location matrix of the inverse-Wishart prior.}
	\item{gamma}{
	A numeric vector of length 1 corresponding to the prior degrees of freedom of the error covariance matrix. The minimum value is m+1, and this is the default value.}
}

\details{
  For technical details of the model, see the accompanying vignette.
}

\value{
The function returns an object of class BVARW, which contains:

  \item{Beta}{
	A matrix of size (m p + 1_c) x m containing the mean posterior estimate of the coefficient matrix, beta.}
	\item{BDraws}{
	An array of size (m p + 1_c) x m x keep which contains the post-burn-in draws of beta.}
	\item{Sigma}{
	A matrix of size m x m containing the mean posterior estimate of the residual covariance matrix, Sigma.}
	\item{SDraws}{
	An array of size m x m x keep which contains the post-burn-in draws of Sigma.}
	\item{IRFs }{
	A four-dimensional object of size irf.periods x m x keep x m containing the impulse response function calculations; the first m refers to the responses to the last m shock. Note: IRF calculations are based on a one standard deviation shock to each variable.}
	\item{data}{
	The data used for estimation.}
	\item{constant}{
	A logical value, TRUE or FALSE, indicating whether the user chose to include a vector of constants in the model.}
}

\author{Keith O'Hara}

\references{
  Koop, Gary and Dimitris Korobilis, ``Bayesian Multivariate Time Series Methods for Empirical Macroeconomics,'' \emph{Mimeo}, 2010.
}

\seealso{
  \code{\link{forecast}}, \code{\link{IRF}}, \code{\link{plot}}.
}

\examples{
\dontrun{
data(BMRVARData)
testbvarw <- BVARW(USMacroData[,2:4],cores=1,c(0.9,0.95,0.95),p=4,constant=T,
                    irf.periods=20,keep=10000,burnin=5000,
                    XiBeta=4,XiSigma=1,gamma=4)
plot(testbvarw,save=F)
IRF(testbvarw,save=F)
forecast(testbvarw,shocks=TRUE,backdata=10,save=FALSE)
}
}

\keyword{BVAR}
\concept{BVAR}