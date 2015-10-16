\name{FREEfit}
\alias{FREEfit}
\alias{FREEfit.formula}
\alias{FREEfit.default}
\title{Function Regression in Ecology and Evolution}
\description{
This is the main call for fitting a model of class \code{FREEfit}. Function \code{FREEfit} has a formula or default method. \code{FREEfit} fits models with function-valued response variables or function-valued predictor variables (but not both). Parameters for these models are represented as continuous functions. Models with function-valued response variables can accommodate clustering variables, which introduce function-valued intercepts for distinct groups within the data (similar to a random intercept in a mixed effects model). The model with a function-valued response variable is:
\deqn{y(m) = \alpha(m) + \Sigma_{i} \beta_{i}(m) x + \Sigma_{q} \gamma_{q, G(i)}(m) + \epsilon(m)},
where \eqn{y(m)} is the response variable as a function of \eqn{m}, \eqn{\alpha(m)} is an intercept function, \eqn{\beta_{i}(m)} is the effect of variable \eqn{i} at each value of \eqn{m}, \eqn{x} is the value of covariate or predictor variable \eqn{i}, \eqn{\gamma_{q, G(i)}(m)} is a function-valued intercept for group \eqn{G(i)} of clustering variable \eqn{q}, where \eqn{G(i)} is the group to which subject \eqn{i} belongs, and \eqn{\epsilon(m)} is the error at each value of \eqn{m}. Models without clustering variables can be fitted using several different methods; models with clustering variables use a custom-built Gibbs sampler.

Models with single-valued (scalar) response variables can accommodate one function-valued predictor variable, one or more clustering variables that introduce group-specific intercepts (similar to random intercepts in a mixed effects model), and one or more scalar predictor variables. The model with a scalar response variable is:
\deqn{y = \alpha + \int \beta(m) x(m) dm + \Sigma_{q} \gamma_{q, G(i)} + \Sigma_{k} \delta_{k} z_{k} + \epsilon},
where \eqn{y} is the response variable, \eqn{\alpha} is the model intercept, \eqn{\beta(m)} is the effect of the function-valued predictor variable, \eqn{x(m)}, on the response variable, \eqn{\gamma_{q, G(i)}} is a random intercept for group \eqn{G(i)} of clustering variable \eqn{q}, where \eqn{G(i)} is the group to which subject \eqn{i} belongs, \eqn{\delta_{k}} is the linear effect of scalar predictor variable, \eqn{z_{k}}, on the response, and \eqn{\epsilon} is the model error. All models with scalar response variables are fitted using a custom-built Gibbs sampler.
}
\usage{
FREEfit(y, ...)
\method{FREEfit}{formula}(formula, data, bins=NULL, groups=NULL, z=NULL, ...)
\method{FREEfit}{default}(y, x, bins = NULL, groups=NULL, z=NULL, ...)
}
\arguments{
  \item{y}{response variable to be passed to \code{FREEfit.method}. Can be function-valued or single-valued}
  \item{formula}{model to be fitted. Of the form y ~ x1 + x2. Interactions are possible, but not tested extensively. Predictor variables can be factor or numeric. If a function-valued predictor is used, then only one predictor variable can be included and the response variable must be single-valued}
  \item{data}{named list containing data for the model. The response variable must be a matrix, with sites in rows, while the predictor variables must be vectors. Predictor variables can be numeric or factors, and are prepared using model.frame()}
  \item{x}{matrix of predictor variables or a matrix of the function-valued predictor variable. Sites must be in rows.}
  \item{bins}{vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{groups}{matrix of clustering variables}
  \item{z}{matrix of scalar predictor variables; only used if the response variable is single-valued}
  \item{\dots}{additional arguments to be passed to \code{\link{FREEspline}} or \code{\link{FREEscalar}}}
}
\details{\code{FREEfit} is the main function for fitting \code{FREEfit} models, and \code{FREEfit.formula} is the main interface.
	
There are many parameters and settings that can be passed to each method; see the appropriate help files for details. These settings require consideration because the different methods will not necessarily perform well for all data under default settings. In many cases the necessary settings might cause slow performance, e.g., many iterations for MCMC methods. See also Yen et al. (2015) for a detailed discussion of performance versus computing time for these methods. We advise running small datasets to test the computing time for each method.
	
The function \code{\link{FREEfitCV}} has been provided for model validation and should be used when determining appropriate settings and for choosing among alternative approaches.}
\value{
An object of class \code{FREEfit} with at least the following components:
  \item{fitted}{fitted y values}
  \item{fitted.sd}{standard deviations of fitted y values}
  \item{observed}{observed y data}
  \item{coefs.mean}{mean values for all fitted coefficients (including the mean function)}
  \item{coefs.sd}{standard deviations of all fitted coefficients (including the mean function)}  
  \item{r2}{r-squared value for the fitted model, based on Pearson R}
  \item{family}{family used for model fitting}
  \item{DIC}{DIC value from the fitted model}
  \item{method}{the method used to fit model}
  \item{call}{If \code{FREEfit.formula} is used: the model call}
  \item{formula}{If \code{FREEfit.formula} is used: the model formula}
}
\references{
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\section{Warning}{BUGS and stan methods will return a warning if any Rhat values are greater than 1.1 (for models where n.chains > 1). This indicates that model convergence might be an issue and should not be ignored. Assessing model convergence and chain mixing is complex and further developments are underway to improve user assessment of MCMC convergence and mixing. Rhat values are returned for each parameter when groups or z are not NULL, along with the log-likelihood at each iteration.}
\seealso{
\code{\link{FREEfit.formula}}, \code{\link{FREEfit.default}},
\code{\link{FREEspline}}, \code{\link{FREEscalar}}, \code{\link{FREEfitCV}}
}
\examples{
# basic model
data <- FREEdataSim(function.class="response", n=10)
data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, n.iters=15)
summary(model)
plot(model)

\dontrun{
# try a model with clustering variables
groups <- sample(1:4, size=nrow(data.list$y), replace=TRUE)
model2 <- FREEfit(y ~ VAR1 + VAR2, data=data.list, groups=groups, n.iters=15)
summary(model2)
plot(model2)

# try a model with a scalar response variable
data <- FREEdataSim(function.class="predictor")
data.list <- with(data, list(y=y.iid, VAR1=x[, , 1]))
z <- rnorm(length(data.list$y))
groups <- sample(1:4, size=length(data.list$y), replace=TRUE)
model3 <- FREEfit(y ~ VAR1, data=data.list, z=z, groups=groups, n.iters=15)
summary(model3)
plot(model3)
}

}
\keyword{function regression}