# This file contains the median similarity imputation method from the
# Goldenberg lab.

medianSimilarity <- function(matrices) {
    # Imputes the missing values in the provided list of matrices
    # with the median values of the respective row.
    #
    # Args:
    #   matrices: A list of matrices that may contain some NAs.
    #
    # Returns:
    #   An imputed version of matrices where all missing values have been replaced
    #   by their appropriate median value.
    m <- length(matrices)
  
    for (v in 1:m) {
        if (is.null(matrices[[v]])) {
            next()
        }
      
        W <- matrices[[v]]
        
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
    
        matrices[[v]] <- W
    }

  
    return(matrices)
}