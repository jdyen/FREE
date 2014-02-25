FREEfitTime <-
function(formula, data=list(), bins=NULL, method=c("fda", "gamboost", "INLA", "stan", "BUGS", "FREE"), n.rep=100, n.cv=10){
  if (length(method) == 1) {
    time.elapsed <- system.time(replicate(n.rep, FREEfitCV(formula=formula, data=data, bins=bins, method=method, n.cv=n.cv, verbose=FALSE)))[["elapsed"]]
  } else {
    time.elapsed <- NULL
    for (i in method) {
      time.temp <- system.time(replicate(n.rep, FREEfitCV(formula=formula, data=data, bins=bins, method=i, n.cv=n.cv, verbose=FALSE)))[["elapsed"]]
      time.elapsed <- c(time.elapsed, time.temp)
      cat("Method: ", i, " is complete.....\n")
      flush.console()
    }
  }
  time.elapsed
}
