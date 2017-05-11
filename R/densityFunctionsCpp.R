# calc_b_splines
#
# Description: updates the B-splines for new values of theta
#
# Arguments:
#
# theta          the locations of the knots for the B-spline basis
# w              the observation points for each value of y
# endpoints      the left- and right-hand endpoints for the B-spline
# degree         the degree of the piecewise polynomial functions for the B-spline

calc_b_splines <- function(theta, w, endpoints=NULL, degree) {
  if (is.numeric(endpoints)) {
    if (length(endpoints) != 2) {
      stop("Error in calc_b_splines: left and right hand end points must be set...",
           call.=FALSE)
    }
  } else {
    if (is.null(endpoints)) {
      endpoints <- c(min(w, na.rm=TRUE) - 1, max(w, na.rm=TRUE) + 1)
    } else {
      stop("Error in calc_b_splines: left and right hand end points must be numeric...",
           call.=FALSE)
    }
  }
  bs.tmp <- bs(x=c(w), knots=theta, degree=degree, Boundary.knots=endpoints)
  b_splines <- array(NA, dim=c(nrow(w), ncol(bs.tmp), ncol(w)))
  for (i in 1:dim(b_splines)[3]) {
    b_splines[, , i] <- bs.tmp[((i - 1) * nrow(w) + 1):(i * nrow(w)), ]
  }
  return(b_splines)
}