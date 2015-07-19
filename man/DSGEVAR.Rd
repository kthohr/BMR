\name{DSGEVAR}

\alias{DSGEVAR}
\alias{DSGEVAR.default}

\title{DSGE-VAR Estimation.}

\description{
  Estimate a DSGE-VAR model.
}

\usage{
DSGEVAR(dsgedata,chains=1,cores=1,lambda=Inf,p=2,
        constant=FALSE,ObserveMat,initialvals,partomats,
        priorform,priorpars,parbounds,parnames=NULL,
        optimMethod="Nelder-Mead",
        optimLower=NULL,optimUpper=NULL,
        optimControl=list(),
        IRFs=TRUE,irf.periods=20,scalepar=1,
        keep=50000,burnin=10000,
        tables=TRUE)
}

\arguments{
  \item{dsgedata}{
	A matrix or data frame of size T x j containing the data series used for estimation. 
  
  Note: in order to identify the structural shocks, there must be the same number of observable series as there are shocks in the DSGE model.}
  \item{chains}{
  A positive integer value indicating the number of MCMC chains to run.}
  \item{cores}{
  A positive integer value indicating the number of CPU cores that should be used for estimation. This number should be less than or equal to the number of chains. 
  
  DO NOT ENTER MORE CORES THAN YOUR COMPUTER CAN SAFELY HANDLE! If in doubt, set cores = 1, which is the default.}
  \item{lambda}{
  The proportion of DSGE dummy data to actual data. Acceptable values lie in the interval j x (p+1)/T to +infinity. (See the vignette for details.)} 
  \item{p}{
  The number of lags to include of each variable. The default value is 2.}
  \item{constant}{
  A logical statement on whether to include a constant vector (intercept) in the model. The default is `FALSE', and the alternative is `TRUE'.}
	\item{ObserveMat}{
	The (m+n+k) x j observable matrix H that maps the state variables to the observable series in the measurement equation.}
	\item{initialvals}{
	Initial values to begin the optimization routine.}
	\item{partomats}{
  This is perhaps the most important function input.
  
  `partomats' should be a function that maps the deep parameters of the DSGE model to the matrices of a solution method, and contains: a k x k matrix labelled `shocks' containing the variances of the structural shocks; a j x 1 matrix labelled `MeasCons' containing any constant terms in the measurement equation; and a j x j matrix labelled `MeasErrs' containing the variances of the measurement errors.}
	\item{priorform}{
	The prior distribution of each parameter.}
	\item{priorpars}{
	The parameters of the prior densities. 
  
  For example, if the user selects a Gaussian prior for a parameter, then the first entry will be the mean and the second its variance.}
	\item{parbounds}{
	The lower- and (where relevant) upper-bounds on the parameter values. `NA' values are permitted.}
	\item{parnames}{
	A character vector containing labels for the parameters.}
	\item{optimMethod}{
  The optimization algorithm used to find the posterior mode. The user may select: the ``Nelder-Mead'' simplex method, which is the default; ``BFGS'', a quasi-Newton method; ``CG'' for a conjugate gradient method; ``L-BFGS-B'', a limited-memory BFGS algorithm with box constraints; or ``SANN'', a simulated-annealing algorithm. 
  
  See \code{\link{optim}} for more details. 
  
	If more than one method is entered, e.g., c(Nelder-Mead, CG), optimization will proceed in a sequential manner, updating the initial values with the result of the previous optimization routine.}
	\item{optimLower}{
	If optimMethod="L-BFGS-B", this is the lower bound for optimization.}
	\item{optimUpper}{
	If optimMethod="L-BFGS-B", this is the upper bound for optimization.}
  \item{optimControl}{
  A control list for optimization. See \code{\link{optim}} for more details.}
	\item{IRFs}{
	Whether to calculate impulse response functions.}
	\item{irf.periods}{
	If IRFs=TRUE, then use this option to set the IRF horizon.}
	\item{scalepar}{
	The scaling parameter, c, for the MCMC run.}
	\item{keep}{
	The number of replications to keep. If keep is set to zero, the function will end with a normal approximation at the posterior mode.}
	\item{burnin}{
	The number of sample burn-in points.}
  \item{tables}{
  Whether to print results of the posterior mode estimation and summary statistics of the MCMC run.}
}

\details{
  For technical details, see the accompanying vignette.
}

\value{
The function returns an object of class DSGEVAR, which contains:
  
  \item{Parameters}{
  A matrix with `keep x chains' number of rows that contains the estimated, post sample burn-in parameter draws.}
  \item{Beta}{
  An array of size (j x p) x m x (keep x chains) which contains the post burn-in draws of beta.}
  \item{Sigma}{
	An array of size j x j x (keep x chains) which contains the post burn-in draws of Sigma.}
	\item{DSGEIRFs}{
	A four-dimensional object of size irf.periods x (m + n + k) x n x (keep x chains) containing the impulse response function calculations for the DSGE model. The first m refers to responses to the last m shock.}
	\item{DSGEVARIRFs}{
	A four-dimensional object of size irf.periods x j x n x (keep x chains) containing the impulse response function calculations for the VAR. The last m refers to the structural shock.}
  \item{parMode}{
  Estimated posterior mode parameter values.}
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
	The data used for estimation.}
}

\author{Keith O'Hara}

\seealso{
  \code{\link{forecast.DSGEVAR}}, \code{\link{IRF.DSGEVAR}}, \code{\link{optim}}, \code{\link{plot.DSGEVAR}}, \code{\link{states.DSGEVAR}}.
}

\examples{
\dontrun{
NKMVAR <- DSGEVAR(dsgedata,chains=1,cores=1,lambda=1,p=4,
                  FALSE,ObserveMat,initialvals,partomats,
                  priorform,priorpars,parbounds,parnames,
                  optimMethod=c("Nelder-Mead","CG"),
                  optimLower=NULL,optimUpper=NULL,
                  optimControl=list(maxit=20000,reltol=(10^(-12))),
                  IRFs=TRUE,irf.periods=5,
                  scalepar=0.28,keep=25000,burnin=25000)
}
}

\keyword{DSGEVAR}
\concept{DSGEVAR}