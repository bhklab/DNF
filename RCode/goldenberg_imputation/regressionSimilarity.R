# This file contains code for performing the regression imputation 
# method from the Goldenberg lab.

regressionSimilarity <- function(matrices) {
    # Imput the missing values for a given matrix in the list of matrices 
    # by creating a linear model from the remaining matrices and using that
    # linear model to predict the missing values.
    # 
    # Args:
    #   matrices: A list of matrices that may contain NAs.
    #
    # Returns:
    #   The provided list of matrices with their missing values imputed.
    m <- length(matrices)
    newMatrices <- vector("list", m)
    
    for (v in 1:m) {
        W <- matrices[[v]]
        
        if (any(is.na(W))) {
            # Replace the diagonal missing values.
            y <- diag(W)
            x <- Reduce(cbind, lapply(matrices[-v], diag))
            x <- as.matrix(x)
            # Replace the missing values in the features with the mean similarity.
            for (k in 1:(m-1)) {
                x[is.na(x[, k]), k] <- mean(x[, k], na.rm=T)
            }
            
            diagonalData <- as.data.frame(cbind(y, x))
            diagonalModel <- lm(y~., diagonalData)
            diag(W)[is.na(y)] <- predict(diagonalModel, diagonalData[is.na(y), ])
            
            # Replace the off-diagonal missing values.
            offDiag <- function(W) {
                diag(W) <- Inf
                W <- W[!is.infinite(W)]
                return(W)
            }
            
            y <- offDiag(W)
            x <- Reduce(cbind, lapply(matrices[-v], offDiag))
            x <- as.matrix(x)
            
            for (k in 1:(m-1)) {
                x[is.na(x[, k]), k] <- mean(x[, k], na.rm=T)
            }
            
            offDiagonalData <- as.data.frame(cbind(y, x))
            offDiagonalModel <- lm(y~., offDiagonalData)
            missingInd <- is.na(W)
            W[missingInd] <- predict(offDiagonalModel, offDiagonalData[is.na(y), ])
        }
        
        newMatrices[[v]] <- W
    }
    
    return(newMatrices)
}