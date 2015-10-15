\name{FREE-package}
\alias{FREE-package}
\alias{FREE}
\docType{package}
\title{Function Regression in Ecology and Evolution}
\description{
\pkg{FREE} fits models where response variables or predictor variables are functions rather than scalar values. These functions draw on several existing packages in R, along with some custom-built MCMC samplers to handle more complex model structures. There is a focus on ecological and evolutionary models, with aids for model interpretation and model validation.
}
\details{
\tabular{ll}{
Package: \tab FREE\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2015-10-14\cr
License: \tab GPL (>= 2)\cr
}
\pkg{FREE} is centred around two main functions: \code{\link{FREEfit}} and \code{\link{FREEfitCV}}. \code{\link{FREEfit}} has a formula interface and takes functional response data (in matrix format, subjects in rows) or single-valued response data (in numeric vector format). Functional response data must be combined with scalar predictor variables (in matrix format, subjects in rows), while single-valued response data can be combined with one functional predictor variable (in matrix format, subjects in rows) and one or more single-valued predictor variables (matrix format, subjects in rows). Models for both functional and single-valued response variables can accommodate clustering variables, which introduce separate intercepts for clusters within the data (similar to random intercepts in a mixed effects model). \code{\link{FREEfit}} fits a model using a method selected by the user (we now support only the \code{FREE} method; see \cite{Yen et al. (2015))}.

\code{\link{FREEfitCV}} performs k-fold cross validation on models, fitted using \code{\link{FREEfit}}. This is useful for model validation and for assessment and tuning of model parameters.
}
\author{
  Jian Yen \cr
  Maintainer: Jian Yen \email{jdl.yen@gmail.com}
}
\references{
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\keyword{ package }
\keyword{ function regression }
\keyword{ ecology }
\keyword{ evolution }
\examples{
# Fit a functional regression model using a Gibbs sampler and plot the results
# using default plot methods for objects of class FREEfit
data <- FREEdataSim(function.class="response", n=10)
data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
groups <- sample(1:4, size=nrow(data.list$y), replace=TRUE)
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, groups=groups, n.iters=15)
summary(model)
plot(model)
}
