medianSimilarity <- function(Wall) {
  m <- length(Wall)
  
  for (v in 1:m) {
    if (is.null(Wall[[v]])) {
        next()
    }
      
    W <- Wall[[v]]
    
    # Calculate the diagonal and off-diagonal medians.
    diagonalMedian <- median(diag(W), na.rm=T)
    offDiagonalMedians <- rep.int(NA, nrow(W))
    for (i in 1:nrow(W)) {
      offDiagonalMedians[i] <- median(W[i, -i], na.rm=T)
    }
    
    # Replace the missing entries with the median entries.
    for (i in 1:nrow(W)) {
      if (is.na(W[i,i])) {
        W[i,i] <- diagonalMedian
      }
      W[i, -i][is.na(W[i, -i])] <- offDiagonalMedians[i]
      W[-i, i][is.na(W[-i, i])] <- offDiagonalMedians[i]
    }
    
    # Replace the remaining missing entries with 0.
    W[is.na(W)] <- 0
    
    Wall[[v]] <- W
  }

  
  return(Wall)
}