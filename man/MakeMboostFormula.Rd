\name{MakeMboostFormula}
\alias{MakeMboostFormula}
\title{MAK FORMULA FOR FUNCTIONAL REGRESSION USING GAMBOOST}
\description{Internal use only: makes a formula for a functional regression model using gamboost. Typically called through \code{\link{FREEfit}} for method gamboost}
\usage{
MakeMboostFormula(n.vars = 2, var.names = NULL, model.int = "bbs", model.pred = "bbs", model.site = "brandom", rand.eff = "bmrf", spatial = FALSE, deg.m.int = 2, df.m.int = 4, diff.m.int = 2, deg.m.pred = 2, df.m.pred = 4, diff.m.pred = 2, df.spat = 6, df.mrf = 100, n.knots = 25)
}
\arguments{
  \item{n.vars}{Number of independent variables}
  \item{var.names}{Vector of variable names}
  \item{model.int}{String describing the model to be used for the intercept function}
  \item{model.pred}{String describing the model to be used for the covariate effects}
  \item{model.site}{String describing the model to be used for the site-based random effects}
  \item{rand.eff}{String to be used for the within-site error structure}
  \item{spatial}{Logical: should random effects be included for spatial location of replicates?}
  \item{deg.m.int}{Degree of the intercept b-spline}
  \item{df.m.int}{Degrees of freedom for the intercept b-spline}
  \item{diff.m.int}{Differences constraint for the intercept b-spline}
  \item{deg.m.pred}{Degree of the b-splines for covariate effects}
  \item{df.m.pred}{Degrees of freedom for the b-splines for covariate effects}
  \item{diff.m.pred}{Differences constraint for the b-splines for covariate effects}
  \item{df.spat}{Degrees of freedom for the bspatial() random effects term}
  \item{df.mrf}{Degrees of freedom for the AR1 error term (implemented as a Markov random field)}
  \item{n.knots}{Specifies the number of knots in the fitted b-splines}
}
\value{
A formula for use by mboost models.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  formula <- MakeMboostFormula()
  formula
}