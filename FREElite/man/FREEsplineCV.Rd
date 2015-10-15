\name{FREEsplineCV}
\alias{FREEsplineCV}
\title{Bayesian function regression using purpose-built MCMC}
\description{\code{FREEsplineCV} runs cross-validation for function regression models with function-valed response variables, one or more scalar predictor variables, and one or more clustering variables. The clustering variables (\code{groups}) define distinct groups that are assigned their own function-valued intercept; this is similar to the use of random intercepts in linear mixed models. The standard method for models with clustering varaibles uses a Gibbs sampler for B-spline-based function regression models, with a first-order autoregressive error structure across bins.}
\usage{
FREEsplineCV(n.cv, y, x, groups, w, degree=3, n_knots_beta=11,
             n_knots_gamma=11, n.iters, n.chains=1, n.thin=1,
             hypers=list(psi_main=0.1, phi_main=0.1, psi_gamma=0.1,
                         phi_gamma=0.1, sigma2_beta=10), ...)
}
\arguments{
  \item{n.cv}{Number of cross-validation folds to use.}
  \item{y}{Numeric vector of response data.}
  \item{x}{Matrix of predictor variable function, sites in rows}
  \item{groups}{Matrix of clustering variables}
  \item{w}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{degree}{The degree of the B spline basis functions used in model fitting}
  \item{n_knots_beta}{Number of knots for B splines for the 'fixed' effects in model fitting}
  \item{n_knots_gamma}{Number of knots for B splines for the 'random' effects in model fitting}
  \item{n.iters}{Total number of MCMC iterations (including burnin)}
  \item{n.chains}{Number of MCMC chains.}
  \item{n.thin}{Thinning rate for MCMC chains; defaults to 1}
  \item{hypers}{Optional: the hyperparameters for the spline parameter variance and for the inverse gamma prior for the overall model variance}
  \item{\dots}{Additional arguments to be passed to \code{FREEspline}}
}
\details{\code{FREEspline} fits a function regression model using a purpose-built Gibbs sampler. The C++ code for this Gibbs sampler is available at the \pkg{FREE} website.}
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
\code{\link{FREEspline}}, \code{\link{FREEfit}}
}
\examples{
\dontrun{
  data <- FREEdataSim(function.class="response", n=10)
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  groups <- sample(1:4, size=nrow(data.list$y), replace=TRUE)
  model <- FREEfitCV(y ~ VAR1 + VAR2, data=data.list, groups=groups,
                     n.iters=15, n.cv=3)
}
}
\keyword{function regression}
\keyword{Bayesian}
\keyword{Gibbs sampler}