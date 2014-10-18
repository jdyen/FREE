\name{FREEinla}
\alias{FREEinla}
\title{Bayesian function regression using INLA}
\description{\code{FREEinla} fits a Bayesian functional regression model for functional responses and one or more scalar covariates. The INLA approach is used, which is a numerical approximation to the Bayesian posterior distribution that avoids the need for MCMC. This method is fast, but relies on several numerical solvers and their associated parameters, which can affect model output.}
\usage{
FREEinla(y, x, bins, family="gaussian", errors="ar1", model.int = "spline",
model.pred = "spline", model.site = "iid", verbose = FALSE, diag.int = 1e-05,
diag.pred = 0.01, prec.prior = c(0.1, 0.001),
control.predictor.set = list(compute = T),
control.compute.set = list(cpo = T, dic = TRUE), ...)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{family}{Family for the observed response variable; all families supported by \code{inla} are supported but only poisson and gaussian errors have been tested}
  \item{errors}{Specify the error structure of within-site residuals; options are "iid" or "ar1"}
  \item{model.int}{String specifying how the intercept function is modelled. spline (default) approximates a spline using a second order random walk. Other options are possible, but have not been tested extensively. See \code{\link{f}}}
  \item{model.pred}{String specifying how the covariate effects are modelled. spline (default) approximates a spline using a second order random walk. Other options are possible, but have not been tested extensively. See \code{\link{f}}}
  \item{model.site}{String specifying how site-based random effects are modelled. iid (default) gives a random intercept to each site. Other options are possible, but have not been tested extensively and are not recommended. See \code{\link{f}}}
  \item{verbose}{Logical: passed directly to \code{\link{inla}}. Should model fitting details be printed directly to the console during model fitting?}
  \item{diag.int}{Constant added to the diagonal of the precision matrix for the intercept function. Can be altered to change smoothness of the splines and to improve model performance.}
  \item{diag.pred}{Constant added to the diagonal of the precision matrix for the covariate effects. Can be altered to change smoothness of the splines and to improve model performance.}
  \item{prec.prior}{Prior for the prior precision for the covariate effects (loggamma parameters). See \code{\link{inla.models}} for details.}
  \item{control.predictor.set}{List specifying settings for control.predictor. See \code{\link{inla}}}
  \item{control.compute.set}{List specifying settings for control.compute. See \code{\link{inla}}}
  \item{\dots}{Additional arguments to be passed to \code{inla}}
}
\details{The INLA approach to model fitting is discussed in detail in the \cite{Rue et al.} references (see below) and at \url{r-inla.org/}. This approach is an alternative to MCMC for Bayesian models and is very effective and relatively quick (minutes rather than hours for complicated models). Several parameters do need to be adjusted for reliable model fitting and this is a work in progress from our end. We are trying to make this process as automatic as possible. Careful use of the model validation tools provided by the package \code{FREE} should avoid any major issues (see, e.g., \code{\link{FREEfitCV}})}
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
}
\references{
  \cite{Yen et al. (in press) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution}
  \cr \cr
  \cite{Rue, H. and Martino, S. and Chopin, N. (2009) Approximate Bayesian Inference for latent Gaussian models using Integrated Nested Laplace Approximations, JRSS-series B (with discussion), vol 71, no 2, pp 319-392.}
  \cr \cr
  \cite{Rue, H and Held, L. (2005) Gaussian Markov Random Fields - Theory and Applications Chapman and Hall}

}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEfda}}, \code{\link{FREEgamboost}}, \code{\link{FREEstan}}, 
\code{\link{FREEbugs}}, \code{\link{FREEbugsJump}}, \code{\link{FREEfit}}
}
\examples{
\dontrun{
  # CAN BE SLOW (MINUTES TO RUN)
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="INLA")
  summary(model)
  plot(model)
}
}
\keyword{function regression}
\keyword{INLA}
\keyword{Bayesian}
\keyword{splines}