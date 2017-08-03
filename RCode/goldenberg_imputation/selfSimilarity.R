selfSimilarity <- function(Wall) {
  if (is.matrix(Wall)) {
    isMatrix <- TRUE
    Wall <- list(Wall)
  }
  
  m <- length(Wall)
  
  for (v in 1:m) {
    W <- Wall[[v]]
    
    # Replace missing diagonal entries with 1.
    diag(W)[is.na(diag(W))] <- 1
    
    # Replace missing off-diagonal entries with 0. 
    W[is.na(W)] <- 0
    
    Wall[[v]] <- W
  }
  
  if (isMatrix) {
    Wall <- Wall[[1]]
  }
  
  return(Wall)
}