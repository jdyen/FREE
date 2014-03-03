\name{FREEbugsJump}
\alias{FREEbugsJump}
\title{Function regression using WinBUGS}
\description{Fit a model where the response data are functions using WinBUGS 1.4. Regression coefficients are splines of arbitrary order and continuity and errors can be IID or AR1.}
\usage{
FREEbugsJump(y, x, bins, family="gaussian", errors="ar1", order=3, cont=order,
bugs.file = NULL, n.chains = 3, n.iters = 5000, n.burnin = n.iters/2,
n.thin = 1, debug = FALSE, bugs.dir = NULL, ...)
}
\arguments{
  \item{y}{Matrix of response data, sites in rows.}
  \item{x}{Matrix of predictor variables, sites in rows}
  \item{bins}{Vector of argument values for response functions (see References). Defaults to 1:ncol(y)}
  \item{family}{Family for the observed response variable; only Gaussian errors are implemented currently}
  \item{errors}{Specify the error structure of within-site residuals; options are "ar1" or "iid"}
  \item{order}{Order of the splines use for coefficients. Defaults to 3.}
  \item{cont}{Degree of continuity for coefficient splines (1 is continuous first derivatives, 2 is continuous second derivatives, and so on). Defaults to the order of the spline.}
  \item{bugs.file}{Optional: String naming WinBUGS file (FILENAME.txt), which must be stored in the working directory (or include path in filename). Defaults to NULL, in which case an appropriate file is automatically created.}
  \item{n.chains}{Number of MCMC chains.}
  \item{n.iters}{Total number of MCMC iterations (including burnin)}
  \item{n.burnin}{Number of iterations to be treated as a burnin}
  \item{n.thin}{Thinning rate for MCMC chains; defaults to 1}
  \item{debug}{Logical: should WinBUGS run in debug mode, which will pause and display diagnostic details on errors or after fitting a model}
  \item{bugs.dir}{Path of directory where WinBUGS14.exe is installed}
  \item{\dots}{Additional arguments to be passed to \code{bugs}}
}
\details{\code{FREEbugsJump} uses WinBUGS 1.4 to fit a regression model where the response data are functions. This method uses reversible-jump MCMC to model coefficients as splines with a flexible number of knots and knot locations. The order and degree of continuity of these splines are set by the user. WinBUGS 1.4 and the jump add-in for reversible-jump MCMC must be installed locally for this function to work.
}
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
\code{\link{FREEfda}}, \code{\link{FREEgamboost}}, \code{\link{FREEinla}}, 
\code{\link{FREEstan}}, \code{\link{FREEbugs}}, \code{\link{FREEfit}}
}
\note{
WinBUGS 1.4 can be very slow. This method is faster than \code{\link{FREEbugs}}, but can still take an hour or more for n=100 functions observed at 20 argument values (under default settings: 5000 iterations and 3 chains). This method can perform very well if time is not an issue. If a spline coefficient method is required, methods \code{\link{FREEinla}} and \code{\link{FREEgamboost}} are quicker alternatives.
}
\section{Warning}{\code{FREEbugsJump} returns a warning if any Rhat value are > 1.1 (if n.chains > 1). This means that there has been an issue with model convergence and should not be ignored. Increasing n.iters might improve convergence.}
\examples{
\dontrun{
  data <- FREEdataSim()
  data.list <- with(data, list(y=y.iid, VAR1=x$VAR1, VAR2=x$VAR2))
  model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="BUGSjump", n.iters=1000)
  summary(model)
  plot(model)
}
}
\keyword{WinBUGS}
\keyword{function regression}
\keyword{Bayesian}
\keyword{splines}