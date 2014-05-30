FREEstan <-
function(y, x, bins, family="gaussian", errors="ar1", n.basis=12, stan.file=NULL, n.chains=3, n.iters=2000, n.burnin=n.iters/2, n.thin=1, verbose=FALSE, refresh=max(n.iters/10, 1), stan.model=NULL, ...){
  stop("stan method not supported in OSX version of FREE...", call.=FALSE)
}