\name{FREEfit}
\alias{FREEfit}
\alias{FREEfit.formula}
\alias{FREEfit.default}
\title{Function Regression in Ecology and Evolution}
\description{
This is the main call for fitting a model of class \code{FREEfit}. Function \code{FREEfit} has a formula or default method. \code{FREEfit} fits a model where the response variable is a function and there are one or more scalar predictor variables. Parameters for the functional model are represented as continuous functions, so the model to be fitted is:
\deqn{y(x) = \alpha(x) + \Sigma_{i} \beta_{i}(x) z + \epsilon(x)}
where \eqn{y(x)} is the response variable as a function of \eqn{x}, \eqn{\alpha(x)} is an intercept function, \eqn{\beta_{i}(x)} is the effect of variable \eqn{i} at each value of \eqn{x}, \eqn{z} is the value of covariate or predictor variable \eqn{i} and \eqn{\epsilon(x)} is the error at each value of \eqn{x}.
}
\usage{
FREEfit(y, ...)
\method{FREEfit}{formula}(formula, data, bins=NULL,
method=c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"),
...)
\method{FREEfit}{default}(y, x, bins = 1:ncol(y),
method = c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"),
bootstrap.cis=FALSE, n.bs=50, ...)
}
\arguments{
  \item{y}{Formula or response variable to be passed to \code{FREEfit.method}}
  \item{formula}{Model to be fitted. Of the form y ~ x1 + x2. Interactions are possible, but not tested extensively. Predictor variables can be factor or numeric.}
  \item{data}{Named list containing data for the model. The response variable must be a matrix, with sites in rows, while the predictor variables must be vectors. Predictor variables can be numeric or factors, and are prepared using model.frame().}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{method}{String defining the method to be used. Seven options are available but only six of these are currently implemented. Defaults to gamboost. Methods BUGSspline and BUGSjump can be slow, as can gamboost, INLA and stan if inappropriate settings are used.}
  \item{bootstrap.cis}{For method=gamboost only: Logical: should a bootstrap be used to calculate confidence intervals for gamboost parameter estimates?}
  \item{n.bs}{If bootstrap.cis=TRUE, the number of bootstraps to perform}
  \item{\dots}{Additional arguments to be passed to \code{FREEmethod}. See \code{\link{FREEfda}}, \code{\link{FREEgamboost}}, \code{\link{FREEinla}}, \code{\link{FREEstan}}, \code{\link{FREEbugs}} and \code{\link{FREEbugsJump}} for details.}
}
\details{\code{FREEfit} is the main function for fitting \code{FREEfit} models, and \code{FREEfit.formula} is the main interface.
	
There are many parameters and settings that can be passed to each method; see the appropriate help files for details. These settings require consideration because the different methods will not necessarily perform well for all data under default settings. In many cases the necessary settings might cause slow performance, e.g., many iterations for MCMC methods. In these cases, appropriate alternatives should be used (\code{\link{FREEfda}} for stan models, \code{\link{FREEinla}} for Bayesian spline models, \code{\link{FREEgamboost}} for general spline models). See also Yen et al. (in preparation) for a detailed discussion of performance versus computing time for these methods. We advise running small datasets to test the computing time for each method; although \code{\link{FREEbugsJump}} might be preferred for its statistical advantages, we often find it to perform extremely slowly (hours to weeks for n > 100).
	
The function \code{\link{FREEfitCV}} has been provided for model validation and should be used when determining appropriate settings and for choosing among alternative approaches.}
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
  \item{method}{the method used to fit model}
  \item{call}{If \code{FREEfit.formula} is used: the model call}
  \item{formula}{If \code{FREEfit.formula} is used: the model formula}
}
\references{
  \cite{Yen et al. (in preparation) Functional regression in ecology and evolution}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\section{Warning}{BUGS and stan methods will return a warning if any Rhat values are greater than 1.1 (for models where n.chains > 1). This indicates that model convergence might be an issue and should not be ignored. Assessing model convergence and chain mixing is complex and further developments are underway to improve user assessment of MCMC convergence and mixing.}
\seealso{
\code{\link{FREEfit.formula}}, \code{\link{FREEfit.default}}, \code{\link{FREEfda}}, 
\code{\link{FREEgamboost}}, \code{\link{FREEinla}}, \code{\link{FREEstan}}, 
\code{\link{FREEbugs}}, \code{\link{FREEbugsJump}}, \code{\link{FREEfree}}, \code{\link{FREEfitCV}}
}
\examples{
data <- FREEdataSim()
data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='fda')
plot(model)
\dontrun{
  model2 <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='INLA')
  plot(model2)
  model3 <- FREEfit(y ~ VAR1 * VAR2, data=data.list, method='stan')
  plot(model3)
  model4 <- FREEfit(y ~ VAR1 * VAR2, data=data.list, method='gamboost')
  plot(model4)
}
}
\keyword{function regression}