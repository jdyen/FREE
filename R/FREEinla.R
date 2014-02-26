FREEinla <-
function(y, x, bins, model.int="rw2",  model.pred="rw2", model.site="iid", model.eij="ar1", family.resp="gaussian", order=NULL, diag.int=1e-5, diag.pred=0.01, prec.prior=c(1e-1,1e-3), control.predictor.set=list(compute=T), control.compute.set=list(cpo=T, dic=TRUE), group.mean=FALSE, n.groups=10, group.vars=FALSE, n.groups.var=10, verbose=FALSE){
  data <- ConvertB2A(y, x, bins)
  y.data <- data$y.vector
  X.data <- data$X.vector
  bin.data <- data$bin.vector
  site.data <- data$sites.vector
  n.vars <- ncol(X.data)
  var.names <- colnames(X.data)
  formula <- MakeInlaFormula(n.vars=n.vars, var.names=var.names, model.int=model.int,
                             model.pred=model.pred, diag.int=diag.int, diag.pred=diag.pred,
                             model.site=model.site, prec.prior=prec.prior, 
                             model.eij=model.eij,
                             order=order, group.mean=group.mean, n.groups=n.groups,
                             group.vars=group.vars, n.groups.var=n.groups.var)
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
  mod.inla <- inla(formula, family=family.resp, data=inla.data.file,
                    control.predictor=control.predictor.set,
                    control.compute=control.compute.set, verbose=verbose)
  fitted <- matrix(mod.inla$summary.fitted$mean, nrow=nrow(y), byrow=TRUE)
  fitted.coefs <- names(mod.inla$summary.random[1:{ncol(x) + 1}])
  coef.vals <- NULL
  coef.vals.sd <- NULL
  for (name.use in fitted.coefs) {
    coef.vals <- cbind(coef.vals, mod.inla$summary.random[[name.use]]$mean)
    coef.vals.sd <- cbind(coef.vals.sd, mod.inla$summary.random[[name.use]]$sd)
    if (family.resp == "poisson") {
      coef.vals <- exp(coef.vals)
      coef.vals.sd <- exp(coef.vals.sd)
    }
  }
  colnames(coef.vals) <- c("INTERCEPT", var.names)
  colnames(coef.vals.sd) <- colnames(coef.vals)
  r <- cor(c(fitted), c(y))
  r2 <- r * r
  xIC <- mod.inla$dic$dic
  return(list(fitted=fitted, observed=y, coefs.mean=t(coef.vals), coefs.sd=t(coef.vals.sd),
              r2=r2, family=family, bins=bins, xIC=xIC))
}
