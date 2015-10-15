dens_func <- function(y, n, from, to, kernel, ...) {
  return(density(y, n=n, from=from, to=to, kernel=kernel, ...)$y)
}