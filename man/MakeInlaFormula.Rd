\name{MakeInlaFormula}
\alias{MakeInlaFormula}
\title{CREATE FORMULA FOR INLA FUNCTIONAL REGRESSION MODEL}
\description{Internal use only: this function creates the formula for the INLA model to be fitted. Typically called by \code{\link{FREEfit}} when method INLA is used.}
\usage{
MakeInlaFormula(n.vars = 2, var.names = NULL, model.int = "rw2", model.pred = "rw2", model.site = "iid", model.eij = "ar1", order = NULL, diag.int = 1e-05, diag.pred = 0.1, prec.prior = c(0.1, 0.001), group.mean = FALSE, n.groups = 10, group.vars = FALSE, n.groups.var = 10)
}
\arguments{
  \item{n.vars}{Number of independent variables in the model}
  \item{var.names}{Vector of variable names}
  \item{model.int}{String defining the model to be used for the intercept function}
  \item{model.pred}{String defining the model to be used for the covariate effects}
  \item{model.site}{String defining the model to be used for site-based random effects}
  \item{model.eij}{String defining the model to be used for within-site errors}
  \item{order}{Order of the arP model if model.eij is set to arP.}
  \item{diag.int}{Constant added to the diagonal of the prior precision matrix for the intercept function}
  \item{diag.pred}{Constant added to the diagonal of the prior precision matrix for the covariate effects}
  \item{prec.prior}{Vector giving parameters for the prior precision function for the covariate effects (loggamma function)}
  \item{group.mean}{Logical: should the intercept function argument values be grouped?}
  \item{n.groups}{If group.mean=TRUE, how many groups should be used?}
  \item{group.vars}{Logical: should the covariate effects argument values be grouped?}
  \item{n.groups.var}{If group.vars=TRUE, how many groups should be used?}
}
\value{
A formula for use in INLA models.
}
\references{
  Yen et al. (in preparation) Functional regression in ecology and evolution
}
\author{Jian Yen <jdl.yen@gmail.com>}
\examples{
  \dontrun{
    formula <- MakeInlaFormula()
    formula
  }
}