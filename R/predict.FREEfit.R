predict.FREEfit <-
function(object, newdata=NULL, ...){
  if (is.null(newdata)) {
    y <- fitted(object)
  } else {
    if (!is.null(object$formula)){
      if (is.null(newdata[all.vars(object$formula)[1]][[1]])) {
        newdata[[length(newdata) + 1]] <- rnorm(length(newdata[[1]]))
        names(newdata)[length(newdata)] <- all.vars(object$formula)[1]
      }
      x <- model.matrix(object$formula, newdata)
    } else {
      x <- cbind(rep(1, nrow(newdata)), newdata)
    }
    y <- x %*% coef(object)$mean
  }
  y
}
