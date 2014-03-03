\name{FREEdataSim}
\alias{FREEdataSim}
\title{Simulate functional response data and predictor variables}
\description{\code{FREEdataSim} simulates a dataset of n observations, particularly for use in testing models of class \code{FREEfit}. Covariates and parameters have default values (see below), but all can be set by the user. Noise is added as IID error (within replicate, within bin), AR1 error and MVN error (correlated across all bins).}
\usage{
FREEdataSim(n = 100, n.vars = 2, bins = NULL, y.mean = NULL, x.val = NULL,
beta.mean = NULL, error.mean = 0, error.sd = 1, cov.mean = NULL, cov.sd = NULL)
}
\arguments{
  \item{n}{Number of replicates}
  \item{n.vars}{Number of independent variables}
  \item{bins}{Either the number of bins (scalar value), or a vector of preferred argument values. Defaults to 1:21.}
  \item{y.mean}{Optional: mean function for simulated data. ncol(y.mean) must equal scalar bins or length bins and must equal the length of beta.mean (if !NULL).}
  \item{x.val}{Optional: covariate values with replicates in rows. Number of rows must equal n and number of columns must equal n.vars. Defaults to NULL, which leads to normally distributed covariate values with mean cov.mean and sd cov.sd}
  \item{beta.mean}{Optional: function for covariate effects, with covariates in rows. ncol(beta.mean) must equal scalar bins or length(bins) and nrow(beta.mean) must be less than or equal to n.vars. Defaults to NULL, which leads to a random Asin(Bx + Cpi) curve. Similarly, if nrow(beta.mean) is less than n.vars the missing rows are filled with random sin() curves.}
  \item{error.mean}{Optional: mean value of errors, which are normally distributed. Defaults to 0 (recommended).}
  \item{error.sd}{Optional: standard deviation of errors, which are normally distributed. Defaults to 1.}
  \item{cov.mean}{Optional: used if x.val is not provided. Sets the mean of covariates, which are normally distributed. Needs to have length equal to n.vars (one mean for each covariate). Defaults to NULL, in which case cov.mean is set to a random integer  between 1 and 10.}
  \item{cov.sd}{Optional: sets the sd of covariates, which are normally distributed (if x.val is not provided). Needs to have length equal to n.vars. Defaults to NULL, in which case cov.sd is set to a random number from c(0.5, 1, 1.5, 2)}
}
\details{This is a very flexible function designed to generate several different random datasets. It has been tested only with settings fairly close to default and may exhibit unusual behaviour away from these settings (e.g., large errors relative to mean values).}
\value{
  \item{y.iid}{Simulated curve with IID errors}
  \item{y.ar1}{Simulated curve with AR1 errors}
  \item{y.mvn}{Simulated curve with multivariate normal errors correlated across all of y}
  \item{x}{Covariate values, either input by user or simulated to follow a normal distrbution}
  \item{mean.real}{Simulated mean function, without noise or covariate effects}
  \item{beta.real}{Simulated covariate effects, one row for each covariate}
  \item{error.iid}{Generated errors for the IID curve}
  \item{error.ar1}{Generated errors for the AR1 curve}
  \item{error.mvn}{Generated errors for the MVN curve}
  \item{bins}{The argument values used}
}
\references{
  \cite{Yen et al. (in preparation) Functional regression in ecology and evolution}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\examples{
data <- FREEdataSim()
y.min <- min(data$y.iid[1, ], data$y.ar1[1, ], data$y.mvn[1, ], data$mean.real)
y.max <- max(data$y.iid[1, ], data$y.ar1[1, ], data$y.mvn[1, ], data$mean.real)
plot(data$y.iid[1, ], type='l', ylim=c(y.min, y.max), xlab="Argument", ylab="Response")
lines(data$y.ar1[1, ], lty=2)
lines(data$y.mvn[1, ], lty=3)
lines(data$mean.real, col='red')
data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, method="fda")
plot(model)
}