MakeMboostFormula <-
function(n.vars=2, var.names=NULL, model.int="bbs", model.pred="bbs", model.site="brandom", rand.eff="bmrf", spatial=FALSE, deg.m.int=2, df.m.int=4, diff.m.int=2, deg.m.pred=2, df.m.pred=4, diff.m.pred=2, df.spat=6,  df.mrf=100, n.knots=25){
  if (is.null(var.names)) {
    var.names <- NULL
    for (i in 1:n.vars) {
      var.names <- c(var.names, paste("VAR", i, sep=""))
    }
  }
  int.temp <- paste("y ~ ", model.int , "(bin.int, degree=", deg.m.int, ", df=", df.m.int,
                    ", differences=", diff.m.int, ", knots=",n.knots - 3,", center = FALSE)", sep="")
  var.temp <- NULL
  for (i in 1:n.vars) {
    var.temp <- paste(var.temp, " + ", model.pred, "(bin.int, by=", var.names[i],
                      ", degree=", deg.m.pred, ", df=", df.m.pred,
                      ", differences=", diff.m.pred, ", knots=",n.knots - 1,", center = TRUE)", sep="")
  }
  if (spatial) {
	if (rand.eff == "bbs") {
      append.temp <- paste(" + bspatial(xcoord, ycoord, df=", df.spat, ") + ",
                           model.site, "(SITE) + bbs(bin.int, by=SITE.fact, center=T)",
                           sep="")
    } else {
      if (rand.eff == "bmrf") {
        append.temp <- paste(" + bspatial(xcoord, ycoord, df=", df.spat,
                             ") + bmrf(bin.mrf, bnd=bin.diag, df=", df.mrf, ")",
                             sep="")
      } else {
        append.temp <- paste(" + bspatial(xcoord, ycoord, df=", df.spat, ") + ",
                             model.site, "(SITE) + brandom(SITE, by=bin.fact)",
                             sep="")
      }
    }
  } else {
    if (rand.eff == "bbs") {
      append.temp <- paste(" + ", model.site, "(SITE) + bbs(bin.int, by=SITE.fact, center=T)",
                           sep="")
    } else {
      if (rand.eff == "bmrf") {
        append.temp <- paste(" + bmrf(bin.mrf, bnd=bin.diag, df=", df.mrf, ")", sep="")
      } else {
        if (rand.eff == "rand.lin") {
          append.temp <- paste(" + ", model.site, "(SITE) + brandom(SITE, by=bin.fact)",
                               sep="")
        } else {
          append.temp <- paste(" + ", model.site, "(SITE)", sep="") 
        }
      }
    }
  }
  formula <- as.formula(paste(int.temp, var.temp, append.temp, sep=""))
  return(formula)
}
