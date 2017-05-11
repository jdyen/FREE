FREEinla <-
function(y, x, bins, family="gaussian", errors="ar1", model.int="spline",  model.pred="spline", model.site="iid", verbose=FALSE, diag.int=1e-5, diag.pred=0.01, prec.prior=c(1e-1,1e-3), control.predictor.set=list(compute=T), control.compute.set=list(cpo=T, dic=TRUE), ...){
  data <- ConvertB2A(y, x, bins)
  y.data <- data$y.vector
  X.data <- data$X.vector
  bin.data <- data$bin.vector
  site.data <- data$sites.vector
  n.vars <- ncol(X.data)
  var.names <- colnames(X.data)
  if (model.int == "spline") {
    model.int <- "rw2"
  } else {
    model.int <- model.int
  }
  if (model.pred == "spline") {
    model.pred <- "rw2"
  } else {
    model.pred <- model.pred
  }
  formula <- MakeInlaFormula(n.vars=n.vars, var.names=var.names, model.int=model.int,
                             model.pred=model.pred, diag.int=diag.int, diag.pred=diag.pred,
                             model.site=model.site, prec.prior=prec.prior, 
                             model.eij=errors,
                             order=order)
  inla.data.file <- list(y=y.data, bin.int=bin.data, SITE=site.data)
  for (i in 1:n.vars) {
    var.temp <- list(X.data[, var.names[i]])
    temp.data <- list(bin.data)
    inla.data.file[length(inla.data.file) + 1] <- temp.data
    names(inla.data.file)[[length(inla.data.file)]] <- paste("bin.int", i, sep="")
    inla.data.file[length(inla.data.file) + 1] <- var.temp
    names(inla.data.file)[[length(inla.data.file)]] <- var.names[i]
  }
  inla.data.file[length(inla.data.file) + 1] <- list(bin.data)
  names(inla.data.file)[[length(inla.data.file)]] <- paste("bin.int", (n.vars + 1), sep="")
  mod.inla <- inla(formula, family=family, data=inla.data.file,
                    control.predictor=control.predictor.set,
                    control.compute=control.compute.set, verbose=verbose, ...)
  fitted <- matrix(mod.inla$summary.fitted$mean, nrow=nrow(y), byrow=TRUE)
  fitted.coefs <- names(mod.inla$summary.random[1:{ncol(x) + 1}])
  coef.vals <- NULL
  coef.vals.sd <- NULL
  for (name.use in fitted.coefs) {
    coef.vals <- cbind(coef.vals, mod.inla$summary.random[[name.use]]$mean)
    coef.vals.sd <- cbind(coef.vals.sd, mod.inla$summary.random[[name.use]]$sd)
    if (family == "poisson") {
      coef.vals <- exp(coef.vals)
      coef.vals.sd <- exp(coef.vals.sd)
    }
  }
  colnames(coef.vals) <- c("INTERCEPT", var.names)
  for (i in 1:n.vars) {
    coef.vals[, {i + 1}] <- coef.vals[, {i + 1}] + mod.inla$summary.fixed[var.names[i], 1]
    coef.vals.sd[, {i + 1}] <- coef.vals.sd[, {i + 1}] + mod.inla$summary.fixed[var.names[i], 2]
  }
  colnames(coef.vals.sd) <- colnames(coef.vals)
  r <- cor(c(fitted), c(y))
  r2 <- r * r
  xIC <- mod.inla$dic$dic
  return(list(fitted=fitted, observed=y, coefs.mean=t(coef.vals), coefs.sd=t(coef.vals.sd),
              r2=r2, family=family, bins=bins, xIC=xIC, formula2=formula, model=mod.inla))
}
