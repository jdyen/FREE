\name{FREEscalar}
\alias{FREEscalar}
\title{Bayesian function regression using purpose-built MCMC}
\description{\code{FREEscalar} implements a purpose-built MCMC approach for Bayesian functional regression, where the predictor variables are given by a function and the response variable is a scalar. This is a work in progress and is in testing stage only. The \code{FREEscalar} approach uses a Gibbs sampler for spline-based function regression models, with a 2nd-order random walk process across bins.}
\usage{
FREEscalar(y, x, bins, n.vars, n.bins, family="gaussian", errors="iid",
n.chains=3, n.iters=10000, n.burnin=n.iters/5, n.thin=n.iters/100,
hypers=list(psi=1, phi=100, sigma2_alpha=10, sigma2_beta=10),
inits=list(vari=10, alpha=NULL, betas=NULL))
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{n.vars}{Number of functional predictor variables in the model}
  \item{n.bins}{Number of bins for the functional predictor variables (must be equal for all variables)}
  \item{family}{Family for the observed response variable; only Gaussian errors are implemented currently}
  \item{errors}{Specify the error structure of within-site residuals; options are "ar1" or "iid"}
  \item{n.chains}{Number of MCMC chains.}
  \item{n.iters}{Total number of MCMC iterations (including burnin)}
  \item{n.burnin}{Number of iterations to be treated as a burnin; defaults to n.iters/5}
  \item{n.thin}{Thinning rate for MCMC chains; defaults to 100}
  \item{hypers}{The hyperparameters for the spline parameter variance and for the inverse gamma prior for the overall model variance}
  \item{inits}{Optional: initial values for the variance and the model parameters}
}
\details{\code{FREEscalar} fits a function regression model using a purpose-built Gibbs sampler. The C code for this Gibbs sampler is available at the \pkg{FREE} website. Extensions to this method are in development.}
\value{
Object of class \code{FREEfit} with the following components:
  \item{fitted}{fitted y values}
  \item{observed}{observed y data}
  \item{coefs.mean}{mean values for all fitted coefficients (including the mean function)}
  \item{coefs.sd}{sd values for all fitted coefficients (including the mean function)}
  \item{r2}{r-squared value for the fitted model, based on Pearson R}
  \item{family}{family used for model fitting}
  \item{bins}{vector of argument values for the response function}
  \item{xIC}{DIC value from the fitted model; not currently available}
}
\references{
  \cite{Yen et al. (in preparation) Function regression in ecology and evolution}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEfda}}, \code{\link{FREEgamboost}}, \code{\link{FREEinla}}, 
\code{\link{FREEstan}}, \code{\link{FREEbugs}}, \code{\link{FREEbugsJump}}, 
\code{\link{FREEfit}}, \code{\link{FREEfree}}
}
\examples{
  data <- FREEdataSim(function.class="predictor")
  data.list <- with(data, list(y=y.iid, VAR1=x[, , 1], VAR2=x[, , 2]))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list)
}
\keyword{function regression}
\keyword{functional predictor}
\keyword{Bayesian}
\keyword{Gibbs sampler}