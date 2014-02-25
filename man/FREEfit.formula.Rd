\name{FREEfit.formula}
\alias{FREEfit.formula}
\title{FUNCTION REGRESSION IN ECOLOGY AND EVOLUTION}
\description{This is the formula interace for fitting models of class \code{FREEfit}. Fits a regression model where the response is a function and there are one or more scalar covariates. Several different methods for model fitting are possible (see Details). These methods all require specific R packages to be installed and, in some cases, require additional software to be installed (see Details).}
\usage{
FREEfit.formula(formula, data = list(), bins = NULL, method = c("fda", "gamboost", "INLA", "stan", "BUGSspline", "BUGSjump", "FREE"), ...)
}
\arguments{
  \item{formula}{Model to be fitted. Of the form y ~ x1 + x2. Interactions are possible, but not tested extensively. Predictor variables can be factor or numeric.}
  \item{data}{Named list containing data for the model. The response variable must be a matrix, with sites in rows, while the predictor variables must be vectors. Predictor variables can be numeric or factors, and are prepared using model.frame().}
  \item{bins}{Optional: vector of the argument values at which the response variable has been observed}
  \item{method}{String specifying which method is to be used for model fitting. Defaults to gamboost. fda is the fastest option, while INLA, stan and BUGSjump/BUGSspline are Bayesian implementations.}
  \item{\dots}{Additional arguments to be passed to FREEfitMethod. See ?FREEfda, ?FREEgamboost, ?FREEinla, ?FREEstan, ?FREEbugsJump, ?FREEbugs and ?FREEfree for details.}
}
\details{FREEfit.formula is the main function for fitting FREEfit models. There are many parameters and settings that can be passed to each method; see the appropriate help files for details. These settings require consideration because the different methods will not necessarily perform well for all data under default settings. In many cases the necessary settings might cause slow performance, e.g., many iterations for MCMC methods. In these cases, appropriate alternatives should be used (fda for stan models, INLA for Bayesian spline models, gamboost for general spline models). See also Yen et al. (in preparation) for a detailed discussion of performance versus computing time for these methods. We advise running small datasets to test the computing time for each method; although BUGS might be preferred for its statistical advantages, we often find it to perform extremely slowly (hours to weeks for n > 100).
	
The function FREEfitCV has been provided for model validation and should be used when determining appropriate settings and for choosing among alternative approaches.}
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
  \item{call}{the model call}
  \item{formula}{the model formula}
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\section{Warning}{BUGS and stan methods will return a warning if any Rhat values are greater than 1.1 (for models where n.chains > 1). This indicates that model convergence might be an issue and should not be ignored. Assessing model convergence and chain mixing is complex and further developments are underway to improve user assessment of MCMC convergence and mixing.}
\seealso{
\code{\link{FREEfit.formula}}
\code{\link{FREEfit.default}}
\code{\link{FREEfda}}
\code{\link{FREEgamboost}}
\code{\link{FREEinla}}
\code{\link{FREEstan}}
\code{\link{FREEbugs}}
\code{\link{FREEbugsJump}}
\code{\link{FREEfree}}
\code{\link{FREEfitCV}}
}
\examples{
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='fda')
  plot(model)
  \dontrun{
    model2 <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method='gamboost')
    plot(model2)
    model3 <- FREEfit(y ~ VAR1 * VAR2, data=data.list, method='fda')
    plot(model3)
  }
}
\keyword{function regression}
\keyword{Bayesian}
\keyword{ecology}
\keyword{evolution}