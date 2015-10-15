MakeInlaFormula <-
function(n.vars=2, var.names=NULL, model.int="rw2", model.pred="rw2", model.site="iid", model.eij="ar1", order=NULL, diag.int=1e-5, diag.pred=0.1, prec.prior=c(0.1,1e-3), group.mean=FALSE, n.groups=10, group.vars=FALSE, n.groups.var=10){
  if (is.null(var.names)) {
    var.names <- NULL
    for (i in 1:n.vars) {
      var.names <- c(var.names, paste("VAR", i, sep=""))
    }
  }
  if (group.mean) {
    int.temp <- paste("y ~ -1 + f(inla.group(bin.int, n=", n.groups, "), model='",
                       model.int, "', diagonal=", diag.int, ", constr=F)", sep="")
  } else {
    int.temp <- paste("y ~ -1 + f(bin.int, model='", model.int, "', diagonal=", diag.int, ",
                       constr=F)", sep="")
  }
  var.temp <- NULL
  var.lin.temp <- NULL
  for (i in 1:n.vars) {
    if (group.vars) {
      bin.int.temp <- paste("inla.group(bin.int", i, ", n=", n.groups.var, ")", sep="")
    } else {
      bin.int.temp <- paste("bin.int", i, sep="")
    }
    var.temp <- paste(var.temp, " + f(", bin.int.temp, ",", var.names[i],
                      ", model='", model.pred, "', hyper=list(prec=list(param=c(",
                      paste(prec.prior[1], prec.prior[2], sep=','),
                      "))), diagonal=", diag.pred, ", constr=T)", sep="")
    var.lin.temp <- paste(var.lin.temp, " + ", var.names[i], sep="")
  }
  if (model.eij == "arp") {
    if(is.null(order)){ stop("User must set order of ARp model...", call.=F)}
    append.temp <- paste(" + f(SITE, model='", model.site, "')", " + f(bin.int",
                         (n.vars + 1), ", model='ar', order=", order,
                         ", replicate=SITE)", sep="")
  } else {
  	if (model.eij == "ar1") {
      append.temp <- paste(" + f(SITE, model='", model.site, "')", " + f(bin.int",
                           (n.vars + 1), ", model='", model.eij, "', replicate=SITE)",
                           sep="")
    } else {
      append.temp <- paste(" + f(SITE, model='", model.site, "')", sep="")
    }
  }
  formula <- as.formula(paste(int.temp, var.lin.temp, var.temp, append.temp, sep=""))
  return(formula)
}
