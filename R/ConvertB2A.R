ConvertB2A <-
function(y.matrix, X.matrix, bin.vector, coord.data=NULL){
  # Converts matrix format functional (binned) data to vector format
  #
  # Args:
  #  y.matrix: Response matrix with one row per site and one column per bin
  #  X.matrix: Predictor variable matrix with one row per site and one column
  #            per variable
  #  bin.vector: Vector of the bin midpoints for data y.matrix
  #
  # Returns:
  #  y.vector: Vector format y
  #  X.vector: Matrix version of X.matrix with repeated entries corresponding
  #            to site IDs
  #  sites.vector: Vector of site IDs
  #  bin.vector: Vector of bin values repeated for each site
  if (nrow(y.matrix) != nrow(X.matrix)) {
    stop("y.matrix and X.matrix must have the same number of rows.....",
         call.=FALSE)
  }
  if (ncol(y.matrix) != length(bin.vector)) {
    stop("y.matrix must have the same number of columns as there are bins.....",
         call.=FALSE)
  }
  n.rows <- nrow(y.matrix)
  n.bins <- length(bin.vector)
  y.vector <- as.vector(t(y.matrix))
  sites.vector <- rep(1:n.rows, each=n.bins)
  X.vector <- X.matrix[sites.vector, ]
  X.vector <- as.matrix(X.vector)
  rownames(X.vector) <- 1:nrow(X.vector)
  colnames(X.vector) <- colnames(X.matrix)
  bin.vector <- rep(bin.vector, times=n.rows)
  if (!is.null(coord.data)) {
    coord.vector <- coord.data[sites.vector, ]
    rownames(coord.vector) <- 1:nrow(coord.vector)
    return(list(y.vector=y.vector, X.vector=X.vector, sites.vector=sites.vector,
                bin.vector=bin.vector, coord.vector=coord.vector))
  } else {
    return(list(y.vector=y.vector, X.vector=X.vector, sites.vector=sites.vector,
                bin.vector=bin.vector))
  }
}
