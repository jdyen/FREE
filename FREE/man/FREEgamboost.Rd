\name{FREEgamboost}
\alias{FREEgamboost}
\title{Function regression using model-based boosting}
\description{Fits a functional regression model for functional response variables using the model-based boosting approach of \pkg{mboost}, in particular, \code{\link{gamboost}}. One or more predictor varialbes is expected. This method is best called using \code{\link{FREEfit}}  and setting method = 'gamboost'.}
\usage{
FREEgamboost(y, x, bins, family="gaussian", errors="ar1", model.int = "spline",
model.pred = "spline", model.site = "iid", deg.int = 2, df.int = 8,
diff.int = 2, deg.beta = 2, df.beta = 8, diff.beta = 2, df.spat = 6,
df.mrf = 100, spatial = FALSE, coord.data = NULL, mstop = 1000,
cvm.set = FALSE, weights = NULL, nu.m = 0.01, trace = FALSE, offset = 0, ...)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{family}{Family for the observed response variable; current options are gaussian, poisson, GammaReg and Nbinomial. See \code{\link{gamboost}} for details.}
  \item{errors}{Specify the error structure of within-site residuals; options are "ar1", "iid", "randlin" and "randspline". "randlin" and "randspline" assign a random linear or spline effect to each site.}
  \item{model.int}{String defining the baselearner for the model intercept. Any baselearner from mboost should work but only spline(baselearner bbs) has been tested at this stage.}
  \item{model.pred}{String defining the baselearner for the effects of predictor variables. Any baselearner from mboost should work but only spline (baselearner bbs) has been tested at this stage.}
  \item{model.site}{String defining the baselearner for the site-based random effects. Any baselearner from mboost should work but only iid (baselearner brandom) has been tested at this stage.}
  \item{deg.int}{Degree of the intercept b-spline}
  \item{df.int}{Degrees of freedom for the intercept b-spline}
  \item{diff.int}{Differences constraint for the intercept b-spline}
  \item{deg.beta}{Degree of the b-splines for covariate effects}
  \item{df.beta}{Degrees of freedom for the b-splines for covariate effects}
  \item{diff.beta}{Differences constraint for the b-splines for covariate effects}
  \item{df.spat}{Degrees of freedom for the bspatial() random effects term}
  \item{df.mrf}{Degrees of freedom for the AR1 error term (implemented as a Markov random field)}
  \item{spatial}{Logical: should random spatial effects be included? If spatial=TRUE these effects are modeled using a b-spline with df=df.spat}
  \item{coord.data}{Optional: matrix of spatial locations with replicates in rows (one column each for latitude/longitude, with column names "LATITUDE" and "LONGITUDE"). Used only if spatial=TRUE.}
  \item{mstop}{The number of initial boosting iterations. Can be set through cross validation by setting cvm.set=TRUE.}
  \item{cvm.set}{Logical: should the optimal number of boosting iterations be estimated using cross validation? If cvm.set=TRUE the optimal mstop will print to the console and all returned results will be set at this value.}
  \item{weights}{Optional: weights for the model response}
  \item{nu.m}{Passed to \code{boost_control}. Specifies the step size or shrinkage parameter.}
  \item{trace}{Logical: trace=TRUE generates printout of status information to the console.}
  \item{offset}{sets model offset}
  \item{\dots}{Additional arguments to be passed to \code{gamboost}}
}
\details{\code{\link{gamboost}} does not generate estimates of uncertainty in parameters by default. A bootstrap method has been implemented to address this and is called directly from \code{\link{FREEfit}}.}
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
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\seealso{
\code{\link{FREEfda}}, \code{\link{FREEinla}}, \code{\link{FREEstan}}, 
\code{\link{FREEbugs}}, \code{\link{FREEbugsJump}}, \code{\link{FREEfit}}
}
\examples{
data <- FREEdataSim(function.class="response", n=10)
data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="gamboost", mstop=250)
summary(model)
plot(model)
}
\keyword{function regression}
\keyword{mboost}
\keyword{splines}