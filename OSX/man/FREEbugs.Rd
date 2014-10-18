\name{FREEbugs}
\alias{FREEbugs}
\title{Function regression using WinBUGS}
\description{NOT SUPPORTED FOR OSX VERSION OF FREE. Fits a basis-function-type function regression model using WinBUGS (see references). Data are modelled over a b-spline basis and the coefficients from this basis are regressed on the predictor variables. These coefficients are assumed independent of one another.}
\usage{
FREEbugs(y, x, bins, family="gaussian", errors="ar1", n.basis=12, bugs.file = NULL,
n.chains = 3, n.iters = 2000, n.burnin = n.iters/2, n.thin = 1, debug = FALSE,
bugs.dir = NULL, ...)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{family}{Family for the observed response variable; only Gaussian errors are implemented currently}
  \item{errors}{Specify the error structure of within-site residuals; options are "ar1" or "iid"}
  \item{n.basis}{Degrees of freedom for b-spline basis}
  \item{bugs.file}{Optional: string naming WinBUGS code file (FILENAME.txt), which must be placed in the working directory (or include path in the filename). Defaults to NULL, which leads to this file being created automatically.}
  \item{n.chains}{Number of MCMC chains.}
  \item{n.iters}{Total number of MCMC iterations (including burnin)}
  \item{n.burnin}{Number of iterations to be treated as a burnin}
  \item{n.thin}{Thinning rate for MCMC chains; defaults to 1}
  \item{debug}{Logical: should WinBUGS run in debug mode, which will pause and display diagnostic details on errors or after fitting a model}
  \item{bugs.dir}{Path of directory where WinBUGS14.exe is installed}
  \item{\dots}{Additional arguments to be passed to \code{bugs}}
}
\details{This method can be VERY SLOW and is not recommended for general use. If a basis-function approach is required we recommend stan or fda methods instead. If using this approach, auto-regressive errors in the observed data are possible but IID errors are recommended to improve computing speed.}
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
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\section{Warning }{
\code{FREEbugs} will return a warning if the fitted model has Rhat values > 1.1 (if n.chains > 1). This indicates issues with convergence and should not be ignored.
}
\seealso{
\code{\link{FREEfda}}, \code{\link{FREEgamboost}}, \code{\link{FREEinla}}, 
\code{\link{FREEstan}}, \code{\link{FREEbugsJump}}, \code{\link{FREEfit}}}
\examples{
\dontrun{
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="BUGSspline", n.chains=2, n.iters=100)
  summary(model)
  plot(model)
}
}
\keyword{WinBUGS}
\keyword{function regression}
\keyword{basis functions}
\keyword{Bayesian}