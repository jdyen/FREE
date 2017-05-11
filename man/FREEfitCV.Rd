\name{FREEfitCV}
\alias{FREEfitCV}
\alias{FREEfitCV.formula}
\alias{FREEfitCV.default}
\title{Cross validation for function regression models}
\description{Main interface for cross validation for models of class \code{FREEfit}. Performs k-fold cross validation and returns an estimate of model predictive r2 as well as observed and predicted values. Plot, print and summary methods are available.}
\usage{
FREEfitCV(y, ...)
\method{FREEfitCV}{formula}(formula, data = list(), bins = NULL, groups=NULL,
                            z=NULL, n.cv = 10, verbose = TRUE, n.iters=100,
                            n.burnin=round(n.iters / 5), n.chains=3, par.run=FALSE, ...)
\method{FREEfitCV}{default}(y, x, bins = NULL, groups=NULL, z=NULL,
                            n.cv = 10, verbose = TRUE, n.iters=100,
                            n.burnin=round(n.iters / 5), n.chains=3, par.run=FALSE, ...)
}
\arguments{
  \item{y}{formula or response variable to be passed to \code{FREEfitCV.method}}
  \item{formula}{model to be fitted. Of the form y ~ x1 + x2. Interactions are possible, but not tested extensively. Predictor variables can be factor or numeric}
  \item{data}{named list containing data for the model. The response variable must be a matrix, with sites in rows, while the predictor variables must be vectors. Predictor variables can be numeric or factors, and are prepared using model.frame()}
  \item{x}{matrix of predictor variables, sites in rows}
  \item{bins}{vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{groups}{matrix of clustering variables}
  \item{z}{matrix of scalar predictor variables}
  \item{n.cv}{number of folds for cross validation. Defaults to 10}
  \item{verbose}{logical: should progress be printed out? If TRUE the percentage complete will be printed to the console after each fold has been held out}
  \item{n.iters}{total number of MCMC iterations (including burnin)}
  \item{n.burnin}{number of iterations to be treated as a burnin; defaults to \code{round(n.iters/5)}}
  \item{n.chains}{number of MCMC chains}
  \item{par.run}{logical: if TRUE and \code{Sys.info()["sysname"] == "Darwin"} the \pkg{parallel} package is used to run multiple chains on different cores}	
  \item{\dots}{additional arguments to be passed to \code{FREEscalar} or \code{FREEspline}}
}
\details{
This is the main interface for cross validation for \code{FREEfit} models; typically the formula interface will be used via \code{FREEfitCV} (see \code{\link{FREEfitCV.formula}} and \code{\link{FREEfitCV}}).

Currently all arguments and data must be supplied directly to \code{FREEfitCV} but a method for passing a fitted \code{FREEfit} model directly to \code{FREEfitCV} is in development.}
\value{
An object of class \code{FREEfitCV}
\item{observed}{Observed data}
\item{predicted}{Predicted data when that particular replicate is held out from model fitting}
\item{cv.r2}{r2 value based on cross validation predictions}
\item{call}{the model call}
\item{formula}{the model formula}
}
\references{
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\section{Warning}{Assessing model convergence and chain mixing is complex and further developments are underway to improve user assessment of MCMC convergence and mixing. Rhat values are returned for each parameter when groups or z are not NULL, along with the log-likelihood at each iteration.}
\seealso{
\code{\link{FREEfitCV.default}}, \code{\link{FREEfitCV.formula}}, \code{\link{FREEfit}}, 
\code{\link{FREEscalar}}, \code{\link{FREEspline}}
}
\examples{
\dontrun{
# basic model with cross validation
data <- FREEdataSim(function.class="response", n=10)
data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
modelCV1 <- FREEfitCV(y ~ VAR1 + VAR2, data=data.list,
                      n.iters=15, n.cv=3)
summary(modelCV1)
plot(modelCV1)

# try a model with clustering variables
groups <- sample(1:4, size=nrow(data.list$y), replace=TRUE)
modelCV2 <- FREEfitCV(y ~ VAR1 + VAR2, data=data.list, groups=groups,
                      n.iters=10)
summary(modelCV2)
plot(modelCV2)
}

# try a model with a scalar response variable
data <- FREEdataSim(function.class="predictor")
data.list <- with(data, list(y=y.iid, VAR1=x[, , 1]))
z <- rnorm(length(data.list$y))
groups <- sample(1:4, size=length(data.list$y), replace=TRUE)
modelCV3 <- FREEfitCV(y ~ VAR1, data=data.list, z=z, groups=groups, n.iters=10, n.cv=3)
plot(modelCV3)
 
}
\keyword{function regression}
\keyword{model validation}