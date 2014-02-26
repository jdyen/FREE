\name{FREEfda}
\alias{FREEfda}
\title{Function regression using \pkg{fda}}
\description{This function uses \pkg{fda} to perform function regression for a functional response variable. Wrapper functions are provided for a formula interface and for k-fold cross validation. Features of this function are already provided by \code{\link{fRegress}} in \pkg{fda}, but \code{FREEfda} provides a convenient interface for ecological and evolutionary analyses and allows for direct comparison with alternative approaches. Plot, print, summary, fitted, predict, residuals, and coef methods are provided.}
\usage{
FREEfda(y, x, bins, y.basis = "bspline", nbasis.y = 12, norder.y = 6, beta.basis = "bspline", nbasis.beta = 12, norder.beta = 6, loglam = NULL)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{y.basis}{String defining the basis functions to use for the response variable. Currently only "bspline" (default) and "fourier" are supported.}
  \item{nbasis.y}{Number of basis functions for response variable}
  \item{norder.y}{Order of splines if y.basis="bspline"}
  \item{beta.basis}{String defining the basis functions to use for the coefficient functions. Currently only "bspline" (default) and "fourier" are supported.}
  \item{nbasis.beta}{Number of basis functions for coefficient functions}
  \item{norder.beta}{Order of splines if beta.basis="bspline"}
  \item{loglam}{Value of penalty to be applied during fitting. See ?fRegress and ?fRegressCV for details.}
}
\details{This function is best called using \code{\link{FREEfit}} or \code{\link{FREEfitCV}} with method="fda". This method is extremely quick but the confidence intervals that are returned result from several assumptions about error structure and are not necessarily as reliable as the Bayesian alternatives (see \code{\link{FREEinla}}, \code{\link{FREEstan}} and \code{\link{FREEbugsJump}}). This method cannot be used for generalized error structures. See references for a comparison of performance among these different methods.}
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
\author{Jian Yen \email{jdl.yen@gmail.com>}}
\examples{
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="fda")
  summary(model)
  plot(model)
}
\keyword{function regression}
\keyword{fda}
\keyword{basis functions}