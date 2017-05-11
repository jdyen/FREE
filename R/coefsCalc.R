# coefs_calc
#
# Description: function to calculate the coefficients of a function regression model
#   based on B-spline coefficients and a B-spline basis`
#
# Arguments
#
# beta          the coefficients for the B-spline basis
# theta         the internal knots for the B-spline basis
# degree        degree of the spline functions in the B-spline basis
# grid          values at which the coefficient functions should be calculated
# endpoints     endpoints for the B-spline basis

coefs_calc <- function(beta, theta, degree, grid=NULL, endpoints=NULL) {
  if (is.null(grid)) {
    grid <- seq(min(sapply(theta, min, na.rm=TRUE)), max(sapply(theta, max, na.rm=TRUE)),
                length=max(sapply(theta, length) + 3))
  }
  if (is.null(endpoints)) {
    endpoints <- c(min(grid) - 1, max(grid) + 1)
  }
  coefs.out <- matrix(NA, nrow=length(beta), ncol=length(grid))
  for (i in seq(along=beta)) {
    coefs.out[i, ] <- beta[[i]] %*% t(bs(grid, knots=theta[[i]], degree=degree,
                                         Boundary.knots=endpoints))
  }
  return(coefs.out)
}