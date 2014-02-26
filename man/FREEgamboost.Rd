\name{FREEgamboost}
\alias{FREEgamboost}
\title{Function regression using model-based boosting}
\description{Fits a functional regression model for functional response variables using the model-based boosting approach of \pkg{mboost}, in particular, \code{\link{gamboost}}. One or more predictor varialbes is expected. This method is best called using \code{\link{FREEfit}}  and setting method = 'gamboost'.}
\usage{
FREEgamboost(y, x, bins, coord.data = NULL, model.int = "bbs", model.pred = "bbs", model.site = "brandom", rand.eff = "bmrf", family = "Gaussian", spatial = FALSE, cvm.set = FALSE, weights = NULL, deg.m.int = 2, df.m.int = 8, diff.m.int = 2, deg.m.pred = 2, df.m.pred = 8, diff.m.pred = 2, df.spat = 6, df.mrf = 100, nu.m = 0.01, mstop = 1000, trace = FALSE, offset = 0)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{coord.data}{Optional: matrix of spatial locations with replicates in rows (one column each for latitude/longitude, with column names "LATITUDE" and "LONGITUDE"). Used only if spatial=TRUE.}
  \item{model.int}{String defining the baselearner for the model intercept. Any baselearner from mboost should work but only bbs has been tested at this stage.}
  \item{model.pred}{String defining the baselearner for the effects of predictor variables. Any baselearner from mboost should work but only bbs has been tested at this stage.}
  \item{model.site}{String defining the baselearner for the site-based random effects. Any baselearner from mboost should work but only brandom has been tested at this stage.}
  \item{rand.eff}{String defining the baselearner for the site-by-argument errors (residual error structure). Three options are currently implemented (bbs, bmrf, rand.lin), and give random b-splines, ar1 errors and random linear effects, respectively, all defined within replicates but across the argument. Any other string will set the site-by-argument errors to be IID.}
  \item{family}{String specifying model family. Currently Poisson, GammaReg and Gaussian are implemented. All other strings will default to a Gaussian model.}
  \item{spatial}{Logical: should random spatial effects be included? If spatial=TRUE these effects are modeled using a b-spline with df=df.spat}
  \item{cvm.set}{Logical: should the optimal number of boosting iterations be estimated using cross validation? If cvm.set=TRUE the optimal mstop will print to the console and all returned results will be set at this value.}
  \item{weights}{Optional: weights for the model response}
  \item{deg.m.int}{Degree of the intercept b-spline}
  \item{df.m.int}{Degrees of freedom for the intercept b-spline}
  \item{diff.m.int}{Differences constraint for the intercept b-spline}
  \item{deg.m.pred}{Degree of the b-splines for covariate effects}
  \item{df.m.pred}{Degrees of freedom for the b-splines for covariate effects}
  \item{diff.m.pred}{Differences constraint for the b-splines for covariate effects}
  \item{df.spat}{Degrees of freedom for the bspatial() random effects term}
  \item{df.mrf}{Degrees of freedom for the AR1 error term (implemented as a Markov random field)}
  \item{nu.m}{Passed to \code{boost_control}. Specifies the step size or shrinkage parameter.}
  \item{mstop}{The number of initial boosting iterations. Can be set through cross validation by setting cvm.set=TRUE.}
  \item{trace}{Logical: trace=TRUE generates printout of status information to the console.}
  \item{offset}{sets model offset}
}
\details{\code{\link{gamboost}} does not generate estimates of uncertainty in parameters by default. We are currently working on a bootstrap method for generating these estimates, but this method performs quite slowly and may not be practicable.}
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
  \cite{Yen et al. (in preparation) Functional regression in ecology and evolution}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEfda}}
\code{\link{FREEinla}}
\code{\link{FREEstan}}
\code{\link{FREEbugs}}
\code{\link{FREEbugsJump}}
\code{\link{FREEfit}}
}
\examples{
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="gamboost", mstop=250)
  summary(model)
  plot(model)
}
\keyword{function regression}