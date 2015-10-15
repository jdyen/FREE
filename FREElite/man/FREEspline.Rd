\name{FREEspline}
\alias{FREEspline}
\title{Bayesian function regression using purpose-built MCMC}
\description{\code{FREEspline} fits a function regression models with function-valed response variables, one or more scalar predictor variables, and one or more clustering variables. The clustering variables (\code{groups}) define distinct groups that are assigned their own function-valued intercept; this is similar to the use of random intercepts in linear mixed models. The standard method for models with clustering varaibles uses a Gibbs sampler for B-spline-based function regression models, with a first-order autoregressive error structure across bins.}
\usage{
FREEspline(y, x, groups, w, degree=3, n_knots_beta=5, n_knots_gamma=5,
           n.iters=1000, n.burnin=round(n.iters / 5), n.thin=1,
           n.chains=3, hypers=list(psi_main=0.1, phi_main=0.1, psi_gamma=0.1,
                                   phi_gamma=0.1, sigma2_beta=10),
           inits=list(beta=NULL, gamma=NULL, rho=0.5, theta1=NULL, theta2=NULL,
                      sigma2=1, sigma2_gamma=1),
           par.run=FALSE)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{groups}{Matrix of clustering variables}
  \item{w}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{degree}{The degree of the B spline basis functions used in model fitting}
  \item{n_knots_beta}{Number of knots for B splines for the 'fixed' effects in model fitting}
  \item{n_knots_gamma}{Number of knots for B splines for the 'random' effects in model fitting}
  \item{n.iters}{Total number of MCMC iterations (including burnin)}
  \item{n.burnin}{Number of iterations to be treated as a burnin; defaults to \code{round(n.iters/5)}}
  \item{n.thin}{Thinning rate for MCMC chains; defaults to 100}
  \item{n.chains}{Number of MCMC chains.}
  \item{hypers}{The hyperparameters for the spline parameter variance and for the inverse gamma prior for the overall model variance}
  \item{inits}{Optional: initial values for the variance and the model parameters}
  \item{par.run}{Logical: if TRUE and \code{Sys.info()["sysname"] == "Darwin"} the \pkg{parallel} package is used to run multiple chains on different cores}
}
\details{\code{FREEspline} fits a function regression model using a purpose-built Gibbs sampler. The C++ code for this Gibbs sampler is available at the \pkg{FREE} website.}
\value{
List with the following components:
  \item{fitted}{Fitted y values}
  \item{fitted.sd}{Standard deviation of fitted y values}
  \item{observed}{Observed y data}
  \item{coefs.mean}{Mean values for all fitted coefficients (including the mean function)}
  \item{coefs.sd}{Standard deviations of all fitted coefficients (including the mean function)}
  \item{rand.coefs.mean}{Mean values for all fitted coefficients for clustering variables}
  \item{rand.coefs.sd}{Standard deviations of all fitted coefficients for clustering variables}
  \item{r2}{R-squared value for the fitted model, based on Pearson R}
  \item{family}{Family used for model fitting}
  \item{w}{Vector of argument values for the response function}
  \item{xIC}{DIC value from the fitted model}
  \item{rhats}{Rhat values for the fitted model; values near one indicate good chain mixing}
  \item{sigma2.mean}{Mean value of the overall model variance}
  \item{sigma2.sd}{Standard deviation of the overall model variance}
  \item{sigma2_gamma.mean}{Mean values of the shared variances for each clustering variable}
  \item{sigma2_gamma.sd}{Standard deviation of the share variances for each clustering variable}
  \item{beta.mean}{Mean values for the B-spline coefficients}
  \item{gamma.mean}{Mean values for the B-spline coefficients for clustering variables}
  \item{rho.mean}{Mean value of the autocorrelation parameter for the overall model errors within subjects}
  \item{rho.sd}{Standard deviation of the autocorrelation parameter}
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
  data <- FREEdataSim(function.class="response", n=10)
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  groups <- sample(1:4, size=nrow(data.list$y), replace=TRUE)
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, groups=groups, n.iters=15)
}
\keyword{function regression}
\keyword{Bayesian}
\keyword{Gibbs sampler}
\keyword{slice sampler}