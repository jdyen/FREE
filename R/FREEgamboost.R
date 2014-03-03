FREEgamboost <-
function(y, x, bins, family="gaussian", errors="ar1", model.int="spline", model.pred="spline", model.site="iid", deg.int=2, df.int=8, diff.int=2, deg.beta=2, df.beta=8, diff.beta=2, df.spat=6, df.mrf=100, spatial=FALSE, coord.data=NULL, mstop=1000, cvm.set=FALSE, weights=NULL, nu.m=0.01, trace=FALSE, offset=0, ...){
  if (spatial & is.null(coord.data)) {
    stop("Coordinates are required if spatial=TRUE...", call.=FALSE)
  }
  n.sites <- nrow(y)
  if (spatial) {
    data <- ConvertB2A(y, x, bins, coord.data)
    coord.data <- data$coord.vector
  } else {
    data <- ConvertB2A(y, x, bins)
  }
  y.data <- data$y.vector
  X.data <- data$X.vector
  bin.data <- data$bin.vector
  site.data <- data$sites.vector
  n.vars <- ncol(X.data)
  var.names <- colnames(X.data)
  bin.fact <- factor(bin.data / max(bin.data))
  SITE.fact <- factor(site.data)
  if (errors == "ar1") {
    bin.diag <- MakeMassiveTridiag(MakeTridiag(c(1, rep(2, length(unique(bin.data)) - 2), 1),
                                 -1, -1), n.times=n.sites)
    bin.mrf <- factor(seq(along=bin.data))
    rownames(bin.diag) <- levels(bin.mrf)
    colnames(bin.diag) <- rownames(bin.diag)
  } else {
    bin.diag <- NULL
    bin.mrf <- NULL
  }
  if (model.site == "iid") {
    model.site <- "brandom"	
  } else {
    model.site <- "brandom"
  }
  if (model.int == "spline") {
    model.int <- "bbs"
  } else {
    model.int <- model.int
  }
  if (model.pred == "spline") {
    model.pred <- "bbs"
  } else {
    model.int <- model.int
  }
  formula <- MakeMboostFormula(n.vars=n.vars, var.names=var.names, model.int=model.int,
                               model.pred=model.pred, model.site=model.site, spatial=spatial,
                               deg.m.int=deg.int, df.m.int=df.int, diff.m.int=diff.int,
                               deg.m.pred=deg.beta, df.m.pred=df.beta,
                               diff.m.pred=diff.beta, df.spat=df.spat, df.mrf=df.mrf,
                               rand.eff=errors, n.knots=length(bins))
  if (spatial) {
    mboost.data.file <- list(y=y.data, bin.int=bin.data, bin.fact=bin.fact, SITE=site.data,
                           SITE.fact=SITE.fact, int=rep(1, length(y.data)),
                           xcoord=coord.data$LONGITUDE, ycoord=coord.data$LATITUDE,
                           bin.diag=bin.diag, bin.mrf=bin.mrf)
  } else {
    mboost.data.file <- list(y=y.data, bin.int=bin.data, bin.fact=bin.fact, SITE=site.data,
                           SITE.fact=SITE.fact, int=rep(1, length(y.data)),
                           bin.diag=bin.diag, bin.mrf=bin.mrf)
  }
  for (i in 1:n.vars) {
    var.temp <- list(X.data[, var.names[i]])
    mboost.data.file[length(mboost.data.file) + 1] <- var.temp
    names(mboost.data.file)[[length(mboost.data.file)]] <- var.names[i]
  }
  if (family == "poisson") {
    mod.mboost <- gamboost(formula=formula, data=mboost.data.file,
                           control=boost_control(mstop=mstop, nu=nu.m, trace=trace),
                           weights=weights, family=Poisson(), offset=offset, ...)
  } else {
    if (family == "GammaReg") {
      mod.mboost <- gamboost(formula=formula, data=mboost.data.file,
                             control=boost_control(mstop=mstop, nu=nu.m, trace=trace),
                             weights=weights, family=GammaReg(), offset=offset, ...)
    } else {
      if (family == "NBinomial") {
        mod.mboost <- gamboost(formula=formula, data=mboost.data.file,
                               control=boost_control(mstop=mstop, nu=nu.m, trace=trace),
                               weights=weights, family=NBinomial(), offset=offset, ...)  
      } else {
        mod.mboost <- gamboost(formula=formula, data=mboost.data.file,
                               control=boost_control(mstop=mstop, nu=nu.m, trace=trace),
                               weights=weights, family=Gaussian(), offset=offset, ...)
      }
    }
  }
  if (cvm.set) {
    cvm <- cvrisk(mod.mboost, folds=cv(model.weights(mod.mboost), type="kfold"))
    mod.mboost[mstop(cvm)]
    cat("mstop set to ", mstop(cvm), " following cross validation...\n\n", sep="")
  }
  vars.to.extract <- variable.names(mod.mboost)
  fitted.vars <- coef(mod.mboost, which=names(vars.to.extract))
  if (length(fitted.vars[[1]]) > 0) {
    coef.vals <- fitted.vars[[1]]
  } else {
    coef.vals <- rep(0, length(bins))
  }
  for (i in seq(along=var.names)) {
  	if (length(fitted.vars[[grep(paste(" ", var.names[i], ",", sep=""),
  	                             names(fitted.vars))]]) > 0) {
      coef.vals <- rbind(coef.vals, fitted.vars[[grep(paste(" ", var.names[i], ",", sep=""),
                                                      names(fitted.vars))]])
    } else {
      coef.vals <- rbind(coef.vals, rep(0, length(fitted.vars[[1]])))
    }
  }
  rownames(coef.vals) <- c("INTERCEPT", var.names)
  fitted <- matrix(fitted(mod.mboost), nrow=nrow(y), byrow=TRUE)
  r <- cor(c(fitted), c(y))
  r2 <- r * r
  xIC <- NULL #AIC(mod.mboost)
  return(list(fitted=fitted, observed=y, coefs.mean=coef.vals, coefs.sd=NULL, r2=r2,
              family=family, bins=bins, xIC=xIC))
}
