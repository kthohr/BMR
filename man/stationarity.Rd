\name{stationarity}

\alias{stationarity}

\title{Testing for Stationarity.}

\description{
  ADF test of a unit root and KPSS test of stationarity.
}

\usage{
stationarity(y,KPSSp=4,ADFp=8,print=TRUE)
}

\arguments{
  \item{y}{
	A matrix or data frame containing the series to be used in testing, and should be of size T x m.}
	\item{KPSSp}{
	The number of lags to include for KPSS test.}
	\item{ADFp}{
	The maximum number of (first-differenced) lags to include in the Augmented Dickey-Fuller (ADF) tests. Lag order is determined by minimising the Bayesian information criterion.}
	\item{print}{
	A logical statement on whether the test results should be printed in the output screen. the default is `TRUE'.}
}

\details{
  Remember, the null hypothesis of the ADF test is NOT the same as the KPSS test; 
  the former is a test with a null of a unit root, and the latter is a test with a null of stationarity.
}

\value{
  Returns a list with
  
  \item{KPSS}{
  KPSS test statistics and critical values.}
	\item{ADF}{
	ADF test statistics and critical values.}
	\item{ADFLag}{
	Number of lags used for the ADF tests.}
}

\author{Keith O'Hara}

\references{
  Hamilton, James D., \emph{Time Series Analysis,} Princeton, New Jersey: Princeton University Press, 1994.
  
  Kwiatkowski, Denis, Peter C.B. Phillips, Peter Schmidt, and Yongcheol Shin, ``Testing the Null Hypothesis of Stationarity Against the Alternative of a Unit Root,'' \emph{Journal of Econometrics,} 1992, \emph{54} (1), 159-178.
}

\examples{
\dontrun{
stationarity(USMacroData[,2:4],4,8)
}
}

\keyword{stationarity}
\concept{stationarity}