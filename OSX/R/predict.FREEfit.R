predict.FREEfit <-
function(object, newdata=NULL, ...){
  if (is.null(newdata)) {
    y <- fitted(object)
  } else {
    if (!is.null(object$formula) & !is.matrix(newdata)){
      if (is.null(newdata$`all.vars(object$formula)[1]`)) {
        newdata[[length(newdata) + 1]] <- rnorm(length(newdata[[1]]))
        names(newdata)[length(newdata)] <- all.vars(object$formula)[1]
      }
      x <- model.matrix(object$formula, newdata)
    } else {
      if (is.matrix(newdata)) {
        x <- cbind(rep(1, nrow(newdata)), newdata)
      } else {
        if (is.vector(newdata)) {
          x <- cbind(rep(1, length(newdata)), newdata)
        } else {
          stop("newdata in predict.FREEfit are not vector or matrix...", call.=FALSE)
        }
      }
    }
    y <- x %*% coef(object)$mean
  }
  y
}
