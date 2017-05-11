FREEcvOutput <-
function(x){
  predicted <- matrix(sapply(x, function(x) t(x$predicted)), ncol=ncol(x[[1]]$predicted),
                      byrow=TRUE)
  observed <- matrix(sapply(x, function(x) t(x$observed)), ncol=ncol(x[[1]]$observed),
                     byrow=TRUE)
  return(list(predicted=predicted, observed=observed))
}
