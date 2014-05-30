\name{FREEstan}
\alias{FREEstan}
\title{Bayesian function regression using MCMC}
\description{NOT SUPPORTED FOR OSX VERSION OF FREE. This fits a Bayesian functional regression model for functional responses and one or more scalar covariates. The stan approach is used, which is an MCMC approach different to the more traditional BUGS approaches, relying on their No-U-Turn Sampler, a variant of Hamiltonian Monte Carlo (see \url{http://mc-stan.org/}). This method is relatively fast, due largely to its quick convergence, but relies on several settings that can affect model output.}
\usage{
FREEstan(y, x, bins, family="gaussian", errors="ar1", n.basis=12,
stan.file = NULL, n.chains = 3, n.iters = 2000, n.burnin = n.iters/2,
n.thin = 1, verbose = FALSE, refresh = max(n.iters/10, 1),
stan.model = NULL, ...)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{family}{Family for the observed response variable; only Gaussian errors are implemented currently}
  \item{errors}{Specify the error structure of the basis function coefficients; options are "iid" and "ar1"}
  \item{n.basis}{The basis function degrees of freedom. Defaults to 12, but any positive integer is possible. Recommended values are in the vicinity of \code{length(bins)/2}.}
  \item{stan.file}{Optional: model code for a stan model. This is created by the function \code{stanCodeDefault} but can be replaced if a non-standard model is required (e.g., multivariate normal errors).}
  \item{n.chains}{Number of MCMC chains}
  \item{n.iters}{Number of MCMC iterations, including burnin}
  \item{n.burnin}{Number of iterations to be treated as a burnin}
  \item{n.thin}{Thinning rate for stored iterations, defaults to 1}
  \item{verbose}{Logical: should the stan model print details of fitting progress to the console? Defaults to FALSE, and even under this setting some information is printed to the console. Generally would only need to be set to TRUE if there are issues with compiling code.}
  \item{refresh}{Progress will be printed to the console every \code{refresh} iterations (even if verbose=FALSE). Default setting \code{n.iters/10} should be appropriate, but setting to zero or any negative integer will suppress printing if preferred.}
  \item{stan.model}{Optional: a previously fitted stan model. Supplying this removes the need for stan to re-compile the c++ code, slightly reducing runtime. Fitted stan models are returned by default by \code{FREEstan} (see Value, below).}
  \item{\dots}{Additional arguments to be passed to \code{stan}}
}
\details{The stan method is an alternative to the Gibbs and Metropolis sampling methods implemented in BUGS approaches (see \url{http://mc-stan.org/} for details). It generally is faster than BUGS approaches because it converges quickly. The method implemented here is a basis-function approach that models the functional responses as b-splines and regresses the b-spline coefficients against the predictor variables.
	
Although faster than BUGS, this method still takes minutes to fit models with 100 replicates and 20 argument values (see examples). Progress is printed to the console so the user can assess running time. See also \cite{Yen et al.} (below) for details of the relative performance of spline-based models and basis-function-based models.

Errors currently are treated as IID within sites, but small changes to the stan.code file will allow other error structures (e.g., AR1 or multivariate normal errors).

Several settings can be changed, and should be assessed using the cross validation methods provided (see \code{\link{FREEfitCV}} for details).

To fit these model the package rstan must be installed. In principle this should work on all systems, but we are unable to run rstan on OSX successfully.
}
\value{
Object of class \code{FREEfit} with the following components:
  \item{fitted}{fitted y values}
  \item{observed}{observed y data}
  \item{coefs.mean}{mean values for all fitted coefficients (including the mean function)}
  \item{coefs.sd}{sd values for all fitted coefficients (including the mean function)}
  \item{r2}{r-squared value for the fitted model, based on Pearson R}
  \item{family}{family used for model fitting}
  \item{bins}{vector of argument values for the response function}
  \item{xIC}{DIC value from the fitted model}
  \item{stan.model}{The fitted stan model, which can be passed to future model fittings to avoid re-compiling c++ code}
}
\references{
  \cite{Yen et al. (in preparation) Functional regression in ecology and evolution}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\section{Warning }{
\code{FREEstan} will return a warning if the fitted model has Rhat values > 1.1 (if n.chains > 1). This indicates issues with convergence and should not be ignored.
}
\seealso{
\code{\link{FREEfda}}, \code{\link{FREEgamboost}}, \code{\link{FREEinla}}, 
\code{\link{FREEbugs}}, \code{\link{FREEbugsJump}}, \code{\link{FREEfit}}
}
\examples{
\dontrun{
  # CAN BE SLOW - PROGRESS IS PRINTED TO CONSOLE
  data <- FREEdataSim()
  data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='stan')
  summary(model)
  plot(model)
}
}
\keyword{function regression}
\keyword{stan}
\keyword{basis functions}
\keyword{Bayesian}