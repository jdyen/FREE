\name{FREEscalar}
\alias{FREEscalar}
\title{Bayesian function regression using purpose-built MCMC}
\description{\code{FREEscalar} fits a function regression models with scalar response variables and one function-valued response variable. The \code{FREEscalar} method uses a Gibbs sampler for B-spline-based function regression models. This method can accommodate clustering varaibles, which define distinct groups that are assigned their own function-valued intercept; this is similar to the use of random intercepts in linear mixed models. This method also can include a matrix of scalar predictor variables.}
\usage{
FREEscalar(y, x, z, groups, bins, degree=3, n_knots=8, n.iters=1000,
           n.burnin=round(n.iters / 5), n.thin=1, n.chains=3,
           hypers=list(phi1=0.1, psi1=0.1, phi2=0.1, psi2=0.1, s2_alpha=10,
                       s2_beta=10, s2_delta=10),
           inits=list(alpha=NULL, beta=NULL, gamma=NULL, delta=NULL,
                      sigma2=1, sigma2_gamma=1),
           par.run=FALSE)
}
\arguments{
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
  \item{hypers}{The hyperparameters for the spline parameter variance and for the inverse gamma prior for the overall model variance}
  \item{inits}{Optional: initial values for the variance and the model parameters}
  \item{par.run}{Logical: if TRUE and \code{Sys.info()["sysname"] == "Darwin"} the \pkg{parallel} package is used to run multiple chains on different cores}
}
\details{\code{FREEscalar} fits a function regression model using a purpose-built Gibbs sampler. The C++ code for this Gibbs sampler is available at the \pkg{FREE} website.}
\value{
List with the following components:
  \item{fitted}{Fitted y values}
  \item{fitted.sd}{Standard deviation of fitted y values}
  \item{observed}{Observed y data}
  \item{coefs.mean}{Mean values for the fitted function-valued model parameter}
  \item{coefs.sd}{Standard deviations of the fitted function-valued model parameter}
  \item{r2}{r-squared value for the fitted model, based on Pearson R}
  \item{family}{family used for model fitting}
  \item{xIC}{DIC value from the fitted model}
  \item{rhats}{Rhat values for the fitted model; values near one indicate good chain mixing}
  \item{sigma2.mean}{Mean value of the overall model variance}
  \item{sigma2.sd}{Standard deviation of the overall model variance}
  \item{sigma2_gamma.mean}{Mean values of the shared variances for each clustering variable}
  \item{sigma2_gamma.sd}{Standard deviation of the share variances for each clustering variable}
  \item{alpha.mean}{Mean of the model grand mean}
  \item{alpha.sd}{Standard deviation of grand mean}
  \item{beta.mean}{Mean values for the B-spline coefficients}
  \item{beta.sd}{Standard deviations of B-spline coefficients}
  \item{gamma.mean}{Mean values for the model intercepts for clustering variables}
  \item{gamma.sd}{Standard deviations of the model intercepts for clustering variables}
  \item{delta.mean}{Mean values for the linear effects of scalar predictor variables}
  \item{delta.sd}{Standard deviations of the fitted effects for scalar predictor variables}
  \item{llik_all}{Log-likelihood values for all iterations (including burn-in) and all chains; useful for assessing model convergence and mixing}
}
\references{
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEfit}}
}
\examples{
  data <- FREEdataSim(function.class="predictor")
  data.list <- with(data, list(y=y.iid, VAR1=x[, , 1]))
  z <- rnorm(length(data.list$y))
  groups <- sample(1:4, size=length(data.list$y), replace=TRUE)
  model <- FREEfit(y ~ VAR1, data=data.list, z=z, groups=groups, n.iters=25)
}
\keyword{function regression}
\keyword{functional predictor}
\keyword{Bayesian}
\keyword{Gibbs sampler}