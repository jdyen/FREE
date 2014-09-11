predict.FREEfit <-
function(object, newdata=NULL, ...){
  if (is.null(newdata)) {
    y <- fitted(object)
  } else {
  	if (is.null(object$method) | (object$method != "scalar")) {
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
    } else {
      if (object$method == "scalar") {
        y <- rep(object$intercept.mean, nrow(newdata))
        for (j in 1:nrow(object$coefs.mean)) {
          y <- y + newdata[, ((j - 1) * ncol(object$coefs.mean) + 1):
                 (j * ncol(object$coefs.mean))] %*% coef(object)$mean[j, ]
        }
      }
    }
  }
  y
}
