GetTargetIndices <- function(integrated, data.bench, num.neighbours, top) {
    predictions <- list()
    
    for (k in num.neighbours) {
        drugs.with.targets <- sort(unique(data.bench$MOLECULE_NAME))
        
        temp <- GetNearestNeighbours(k, integrated, drugs.with.targets)
        
        neighbours <- temp$neighbours
        weights <- temp$weights
        
        counts <- CountTargetsFromNeighbours(neighbours, weights, data.bench, 
                                             scale.distance=F, extra.scaling=100)
        top.predictions <- GetTopPredictions(counts, top)    
        predictions[[k]] <- top.predictions
    }
    
    return(predictions)    
}

DetermineTargetDistributions <- function(predictions, data.bench, num.neighbours, top) {
    distributions <- list()
    
    for (k in num.neighbours) {
        distributions[[k]] <- numeric(top)
        temp <- predictions[[k]]
        
        for (drug in names(temp)) {
            pred <- names(temp[[drug]])
            possible.targets <- data.bench$TARGET_NAME[data.bench$MOLECULE_NAME == drug]
            
            for (j in 1:top) {
                if (pred[j] %in% possible.targets) {
                    distributions[[k]][j] <- distributions[[k]][j] + 1
                    break
                }
            }
        }   
    }
    
    return(distributions)
}

CreateTargetDistDataForGGPlot <- function(distributions, accuracies, num.neighbours, top) {
    for (k in num.neighbours) {
        if (k == num.neighbours[1]) {
            df.data <- cbind(rep(k, top), 1:top, distributions[[k]])
        } else {
            temp <- cbind(rep(k, top), 1:top, distributions[[k]])
            df.data <- rbind(df.data, temp)
        }
    }
    
    colnames(df.data) <- c("k", "index", "counts")
    df.data <- as.data.frame(df.data)
    temp <- paste("K:", df.data$k, "\n acc:", accuracies, sep = "")
    df.data$k <- factor(temp, levels=(temp[!duplicated(temp)]))
    df.data$index <- as.factor(df.data$index)
    
    df.data <- cbind(df.data, paste(df.data$k, accuracies, sep=" - "))
    colnames(df.data)[4] <- "acc"
    
    return(df.data)
}