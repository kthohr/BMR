\name{DSGEVAR}

\alias{DSGEVAR}
\alias{plot.DSGEVAR}

\title{DSGE-VAR Estimation.}

\description{
  Estimating a DSGE-VAR model.
}

\usage{
DSGEVAR(dsgedata,lambda=Inf,p=2,ObserveMat,initialvals,partomats,priorform,
        priorpars,parbounds,parnames=NULL,
        optimMethod="Nelder-Mead",optimLower=NULL,optimUpper=NULL,
        optimControl=list(),
        DSGEIRFs=TRUE,irf.periods=20,
        scalepar=1,keep=50000,burnin=10000)
\method{plot}{DSGEVAR}(obj,BinDenom=40,save=FALSE,height=13,width=13)
}

\arguments{
  \item{dsgedata}{
	A matrix or data frame of size T x j containing the relevant series.} 
  \item{lambda}{
  The proportion of DSGE dummy data to actual data. See the vignette for details.} 
  \item{p}{
  The number of lags to include of each variable. The default value is 2.}
	\item{ObserveMat}{
	The (m+n+k) x j observable matrix H that maps the state variables to the to observed series in the measurement equation.}
	\item{initialvals}{
	The initial values to begin the optimisation routine.}
	\item{partomats}{
	This is possibly the most important component the user supplies to the function. It should be a function that maps the deep parameters to the matrices of Uhlig's solution, A through N, and an additional matrix labelled `shocks', which should be of size k x k.}
	\item{priorform}{
	The prior distribution of each parameter.}
	\item{priorpars}{
	The parameters of the relevant prior densities. For example, if the user selects a Gaussian distribution for a parameter, then the first entry will be the mean and the second being the variance.}
	\item{parbounds}{
	The lower- and, where relevant, upper-bounds on the parameter values.}
	\item{parnames}{
	Parameter names. This is useful to add for graphing densities later.}
	\item{optimMethod}{
	The optimisation algorithm used to find the posterior mode. The user may select: the "Nelder-Mead" simplex method, which is the default; "BFGS", a quasi-Newton method; "CG" for a conjugate gradient method; "L-BFGS-B", a limited-memory BFGS algorithm with box constraints; or "SANN", a simulated-annealing algorithm. See ?optim for more details.}
	\item{optimLower}{
	If optimMethod="L-BFGS-B", this is the lower bound for optimisation.}
	\item{optimUpper}{
	If optimMethod="L-BFGS-B", this is the upper bound for optimisation.}
  \item{optimControl}{
  A control list for optimisation. See ?optim in R for more details.}
	\item{IRFs}{
	Whether to calculate impulse response functions.}
	\item{irf.periods}{
	If IRFs=TRUE, then this option sets the horizon of the IRFs.}
	\item{scalepar}{
	The scaling parameter, c, for the MCMC run.}
	\item{keep}{
	The number of replications to keep.}
	\item{burnin}{
	The number of sample burn-in points.}
}

\details{
  For technical details, see the accompanying vignette.
}

\value{
The function returns an object of class DSGEVAR, which contains:

  \item{Parameters}{
	A matrix with `keep' number of rows that contain the estimated, post sample burn-in parameter draws.}
	\item{IRFs}{
	The IRFs (if any), based on the posterior parameter draws.}
	\item{ModeParamTrans}{
	The values of the transformed parameters at the posterior mode.}
	\item{ModeHessian}{
	The Hessian computed at the posterior mode.}
	\item{AcceptanceRate}{
	The acceptance rate of the chain.}
	\item{ObserveMat}{
	The user-supplied H matrix from the Kalman filter recursion.}
	\item{data}{
	The data used in estimation.}
}

\author{Keith O'Hara}

\seealso{
  \code{\link{plot}}, \code{\link{optim}}
}

\examples{
\dontrun{
#
NKMVAR <- DSGEVAR(dsgedata,1,2,ObserveMat,initialvals,partomats,
                  priorform,priorpars,parbounds,parnames,
                  optimMethod="CG",optimLower=NULL,optimControl=list(maxit=20000,reltol=(10^(-12))),
                  optimUpper=NULL,IRFs=TRUE,irf.periods=5,
                  scalepar=0.4,keep=10000,burnin=10000)
}
}

\keyword{DSGEVAR}
\concept{DSGEVAR}