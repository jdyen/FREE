\name{FREEdataSim}
\alias{FREEdataSim}
\title{Simulate functional response or predictor data and all other associated variables}
\description{\code{FREEdataSim} simulates a dataset of n observations, particularly for use in testing models of class \code{FREEfit}. Either functional responses (default) or functional predictors can be simulated. Covariates and parameters have default values (see below), but all can be set by the user. For functional response models, noise is added as IID error (within replicate, within bin), AR1 error and MVN error (correlated across all bins). For functional predictor models, all noise follows a Gaussian IID process.}
\usage{
FREEdataSim(function.class = c("response", "predictor"), n = 100, n.vars = 2,
            bins = NULL, y.mean = NULL, x.val = NULL, beta.mean = NULL,
            error.mean = 0, error.sd = 1, cov.mean = NULL, cov.sd = NULL,
            struc.err.mean = 0, struc.err.sd = 1)
}
\arguments{
  \item{function.class}{should the simulated data have functional responses or functional predictors? Defaults to "response", which generates data with a functional response variable}
  \item{n}{number of replicates}
  \item{n.vars}{number of independent variables}
  \item{bins}{either the number of bins (scalar value), or a vector of preferred argument values. Defaults to 1:21}
  \item{y.mean}{optional: mean function for simulated data. ncol(y.mean) must equal scalar bins or length bins and must equal the length of beta.mean (if !NULL)}
  \item{x.val}{optional: covariate values with replicates in rows. Number of rows must equal n and number of columns must equal n.vars. Defaults to NULL, which leads to normally distributed covariate values with mean cov.mean and sd cov.sd}
  \item{beta.mean}{optional: function for covariate effects, with covariates in rows. ncol(beta.mean) must equal scalar bins or length(bins) and nrow(beta.mean) must be less than or equal to n.vars. Defaults to NULL, which leads to a random Asin(Bx + Cpi) curve. Similarly, if nrow(beta.mean) is less than n.vars the missing rows are filled with random sin() curves}
  \item{error.mean}{optional: mean value of errors, which are normally distributed. Defaults to 0 (recommended)}
  \item{error.sd}{optional: standard deviation of errors, which are normally distributed. Defaults to 1}
  \item{cov.mean}{optional: used if x.val is not provided. Sets the mean of covariates, which are normally distributed. Needs to have length equal to n.vars (one mean for each covariate). Defaults to NULL, in which case cov.mean is set to a random integer  between 1 and 10}
  \item{cov.sd}{optional: sets the sd of covariates, which are normally distributed (if x.val is not provided). Needs to have length equal to n.vars. Defaults to NULL, in which case cov.sd is set to a random number from c(0.5, 1, 1.5, 2)}
  \item{struc.err.mean}{optional: mean of structural (function-wide) errors, which are normally distributed. Defaults to 0}
  \item{struc.err.sd}{optional: standard deviation of structural (function-wide) errors, which are normally distributed. Defaults to 1}
}
\details{This is a very flexible function designed to generate several different random datasets. It has been tested only with settings fairly close to default and may exhibit unusual behaviour away from these settings (e.g., large errors relative to mean values).}
\value{
  \item{y.iid}{simulated curve with IID errors}
  \item{y.ar1}{simulated curve with AR1 errors}
  \item{y.mvn}{simulated curve with multivariate normal errors correlated across all of y}
  \item{x}{covariate values, either input by user or simulated to follow a normal distrbution}
  \item{mean.real}{simulated mean function, without noise or covariate effects}
  \item{beta.real}{simulated covariate effects, one row for each covariate}
  \item{error.iid}{generated errors for the IID curve}
  \item{error.ar1}{generated errors for the AR1 curve}
  \item{error.mvn}{generated errors for the MVN curve}
  \item{bins}{the argument values used}
}
\references{
  \cite{Yen et al. (2015) Function regression in ecology and evolution: FREE. Methods in Ecology and Evolution, 6: 17-26.}
}
\author{Jian Yen \email{jdl.yen@gmail.com}}
\examples{
data <- FREEdataSim(function.class="response", n=15)
y.min <- min(data$y.iid[1, ], data$y.ar1[1, ], data$y.mvn[1, ], data$mean.real)
y.max <- max(data$y.iid[1, ], data$y.ar1[1, ], data$y.mvn[1, ], data$mean.real)
plot(data$y.iid[1, ], type='l', ylim=c(y.min, y.max), xlab="Argument", ylab="Response")
lines(data$y.ar1[1, ], lty=2)
lines(data$y.mvn[1, ], lty=3)
lines(data$mean.real, col='red')
data.list <- list(y=data$y.iid, VAR1=data$x$VAR1, VAR2=data$x$VAR2)
model <- FREEfit(y ~ VAR1 + VAR2, data=data.list, n.iters=15)
plot(model)
}