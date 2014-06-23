\name{EDSGE}

\alias{EDSGE}
\alias{plot.EDSGE}

\title{DSGE Estimation.}

\description{
  Estimating a DSGE model.
}

\usage{
EDSGE(dsgedata,chains=1,cores=1,
      ObserveMat,initialvals,partomats,
      priorform,priorpars,parbounds,parnames=NULL,
      optimMethod="Nelder-Mead",
      optimLower=NULL,optimUpper=NULL,
      optimControl=list(),
      DSGEIRFs=TRUE,irf.periods=20,
      scalepar=1,keep=50000,burnin=10000)
\method{plot}{EDSGE}(obj,BinDenom=40,save=FALSE,height=13,width=13)
}

\arguments{
  \item{dsgedata}{
	A matrix or data frame of size T x j containing the relevant series.}
  \item{chains}{
  A positive integer value indicating the number of MCMC chains to run.}
  \item{cores}{
  A positive integer value indicating the number of CPU cores that should be used in estimation. This number should be less than or equal to the number of chains. DO NOT ENTER MORE CORES THAN YOUR COMPUTER CAN HANDLE! If in doubt, set cores = 1, which is the default.}
	\item{ObserveMat}{
	The (m+n+k) x j observable matrix H that maps the state variables to the to observed series in the measurement equation.}
	\item{initialvals}{
	The initial values to begin the optimisation routine.}
	\item{partomats}{
	This is perhaps the most important input to the function. It should be a function that maps the deep parameters to the matrices of Uhlig's solution, A through N, and an additional matrix labelled `shocks', which should be of size k x k and contain the variances of the structural shocks.}
	\item{priorform}{
	The prior distribution of each parameter.}
	\item{priorpars}{
	The parameters of the relevant prior densities. For example, if the user selects a Gaussian distribution for a parameter, then the first entry will be the mean and the second being the variance.}
	\item{parbounds}{
	The lower- and, where relevant, upper-bounds on the parameter values.}
	\item{parnames}{
	Parameter names. This is useful to add for graphing densities later.}
	\item{optimMethod}{
	The optimisation algorithm used to find the posterior mode. The user may select: the ``Nelder-Mead'' simplex method, which is the default; ``BFGS'', a quasi-Newton method; ``CG'' for a conjugate gradient method; ``L-BFGS-B'', a limited-memory BFGS algorithm with box constraints; or ``SANN'', a simulated-annealing algorithm. See ?optim in R for more details. 
  
	If more than one method is entered, e.g., (``Nelder-Mead'', ``CG''), the function will run optimization in a sequential manner, updating the initial values with the result of the previous optimisation routine.}
	\item{optimLower}{
	If optimMethod="L-BFGS-B", this is the lower bound for optimisation.}
	\item{optimUpper}{
	If optimMethod="L-BFGS-B", this is the upper bound for optimisation.}
  \item{optimControl}{
  A control list for optimisation. See ?optim in R for more details.}
	\item{DSGEIRFs}{
	Whether to calculate impulse response functions.}
	\item{irf.periods}{
	If DSGEIRFs=TRUE, then this option sets the horizon of the IRFs.}
	\item{scalepar}{
	The scaling parameter, c, for the MCMC run.}
	\item{keep}{
	The number of replications to keep. If keep is set to zero, the function will end with the normal approximation at the posterior mode.}
	\item{burnin}{
	The number of sample burn-in points.}
}

\details{
  For technical details, see the accompanying vignette.
}

\value{
The function returns an object of class EDSGE, which contains:

  \item{Parameters}{
	A matrix with `keep x chains' number of rows that contain the estimated, post sample burn-in parameter draws.}
  \item{parMode}{
  The values of the parameters at the posterior mode.}
	\item{ModeHessian}{
	The Hessian computed at the posterior mode for the transformed parameters.}
  \item{logMargLikelihood}{
  The log marginal likelihood from a Laplacian approximation at the posterior mode.}
	\item{IRFs}{
	The IRFs (if any), based on the posterior parameter draws.}
	\item{AcceptanceRate}{
	The acceptance rate of the chain(s).}
  \item{RootRConvStats}{
  Gelman's sqrt{R}-between-chain convergence statistics for each parameter. A value close 1 would signal convergence.}
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
ObserveMat<-rbind(0,0,1,0,0,0,0)
#
initialvals <- c(0.28,0.015,1,0.90,1)
parnames <- c("Alpha","Delta","Eta","Rho","SigmaA")
priorform <- c("Normal","Normal","Normal","Beta","IGamma")
priorpars <- t(rbind(c(0.30, 0.015,   1,   4,  2),
                     c(0.05^2, 0.002^2, 0.2^2, 1.5,  1)))
parbounds <- t(rbind(c(NA, NA, NA,   0.5,   0),
                     c(NA, NA, NA, 0.999,  NA)))
#
RBCest <- EDSGE(dsgedata,1,1,ObserveMat,initialvals,partomats,
           priorform,priorpars,parbounds,parnames,
           optimMethod="BFGS",optimLower=NULL,optimUpper=NULL,
           DSGEIRFs=TRUE,irf.periods=30,scalepar=0.75,keep=50000,burnin=20000)
}
}

\keyword{EDSGE}
\concept{EDSGE}