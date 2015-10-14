\name{FREEscalarCV}
\alias{FREEscalarCV}
\title{Bayesian function regression using purpose-built MCMC}
\description{\code{FREEscalarCV} runs cross-validation for function regression models with scalar response variables and one function-valued response variable. The \code{FREEscalar} method uses a Gibbs sampler for B-spline-based function regression models. This method can accommodate clustering varaibles, which define distinct groups that are assigned their own function-valued intercept; this is similar to the use of random intercepts in linear mixed models. This method also can include a matrix of scalar predictor variables.}
\usage{
FREEscalarCV(n.cv=10, y, x, z, groups, bins, 
             degree=3, n_knots=8, n.iters=1000,
             n.burnin=round(n.iters / 5), n.thin=1, n.chains=3,
             hypers=list(phi1=0.1, psi1=0.1, phi2=0.1, psi2=0.1, s2_alpha=10,
                         s2_beta=10, s2_delta=10),
             inits=list(alpha=NULL, beta=NULL, gamma=NULL, delta=NULL,
                        sigma2=1, sigma2_gamma=1), ...)
}
\arguments{
  \item{n.cv}{Number of cross-validation folds to use.}
  \item{y}{Numeric vector of response data.}
  \item{x}{Matrix of predictor variable function, sites in rows}
  \item{z}{Matrix of scalar predictor variables}
  \item{groups}{Matrix of clustering variables}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{degree}{The degree of the B spline basis functions used in model fitting}
  \item{n_knots}{Number of knots for B splines for the function-valued effects in model fitting}
  \item{n.iters}{Total number of MCMC iterations (including burnin)}
  \item{n.burnin}{Number of iterations to be treated as a burnin; defaults to \code{round(n.iters/5)}}
  \item{n.thin}{Thinning rate for MCMC chains; defaults to 1}
  \item{n.chains}{Number of MCMC chains.}
  \item{hypers}{Optional: the hyperparameters for the spline parameter variance and for the inverse gamma prior for the overall model variance}
  \item{inits}{Optional: initial values for the variance and the model parameters}
  \item{\dots}{Additional arguments to be passed to \code{FREEscalar}}
}
\details{\code{FREEscalar} fits a function regression model using a purpose-built Gibbs sampler. The C++ code for this Gibbs sampler is available at the \pkg{FREE} website.}
\value{
Object of class \code{FREEfitCV} with the following components:
  \item{predicted}{predicted y values}
  \item{observed}{observed y data}
}
\references{
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEscalar}}, \code{\link{FREEfit}}
}
\examples{
  data <- FREEdataSim(function.class="predictor")
  data.list <- with(data, list(y=y.iid, VAR1=x[, , 1]))
  z <- rnorm(length(data.list$y))
  groups <- sample(1:4, size=length(data.list$y), replace=TRUE)
  model <- FREEfitCV(y ~ VAR1, data=data.list, z=z, groups=groups, n.iters=50)
}
\keyword{function regression}
\keyword{Bayesian}
\keyword{Gibbs sampler}