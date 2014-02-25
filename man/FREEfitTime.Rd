\name{FREEfitTime}
\alias{FREEfitTime}
\title{COMPUTING TIME FOR FREEfit METHODS}
\description{Runs FREEfitCV for one or more FREEfit methods and returns the computational time for each method.}
\usage{
FREEfitTime(formula, data = list(), bins = NULL, method = c("fda", "gamboost", "INLA", "stan", "BUGS", "FREE"), n.rep = 100, n.cv = 10)
}
\arguments{
  \item{formula}{Model to be fitted. Of the form y ~ x1 + x2. Interactions are possible, but not tested extensively. Predictor variables can be factor or numeric.}
  \item{data}{Named list containing data for the model. The response variable must be a matrix, with sites in rows, while the predictor variables must be vectors. Predictor variables can be numeric or factors, and are prepared using model.frame().}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{method}{String defining the method to be used. Seven options are available but only six of these are currently implemented. Defaults to gamboost. Methods BUGSspline and BUGSjump can be slow, as can gamboost, INLA and stan if inappropriate settings are used.}
  \item{n.cv}{Number of folds for cross validation. Defaults to 10.}
  \item{n.rep}{Number of times to replicate the FREEfitCV computation}
}
\details{
FOR DIAGNOSTIC USE ONLY - NOT FOR GENERAL USE. CAN BE VERY SLOW SO RUN AT OWN RISK.
}
\value{
\item{time.elapsed}{Scalar or vector reporting the time taken for each different method}
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  \dontrun{
    data <- FREEdataSim()
    data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
    modelTime <- FREEfitTime(y ~ VAR1 + VAR2, data=data.list, method='fda', n.rep=10, n.cv=5)
    modelTime
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
