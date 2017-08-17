GetTargetPerformance <- function(top.predictions, data.bench) {
    all.target.names <- sort(unique(data.bench$TARGET_NAME))
    
    target.performance <- numeric(length(all.target.names))
    names(target.performance) <- all.target.names
    
    for (drug in names(top.predictions)) {
        predictions <- names(top.predictions[[drug]])
        
        true.targets <- data.bench$TARGET_NAME[data.bench$MOLECULE_NAME == drug]
        
        if (length(intersect(predictions, true.targets)) > 0) {
            target.performance[true.targets] <- target.performance[true.targets] + 1
        } 
    }

    return(target.performance)    
}

GetMisclassifiedDrugsPerTarget <- function(top.predictions, data.bench) {
    target.misclassifications <- list()
    
    for (drug in names(top.predictions)) {
        predictions <- names(top.predictions[[drug]])
        
        true.targets <- data.bench$TARGET_NAME[data.bench$MOLECULE_NAME == drug]
        
        if (length(intersect(predictions, true.targets)) == 0) {
            for (target in true.targets) {
                if (is.null(target.misclassifications[[target]])) {
                    target.misclassifications[[target]] <- list()
                }
                
                target.misclassifications[[target]][[drug]] <- drug
            }
        }
    }
    
    return(target.misclassifications)
}

CreateDataForBadTargetsPlot <- function(bad.targets) {
    bad.targets.counts <- unlist(lapply(bad.targets, length))
    
    df.data <- data.frame(network=names(bad.targets.counts), 
                          counts=bad.targets.counts, 
                          stringsAsFactors = F)
    
    return(df.data)
}

MergeTargets <- function(bad.targets, new.targets) {
    for (target in names(new.targets)) {
        if (target %in% names(bad.targets)) {
            bad.targets[target] <- bad.targets[target] + new.targets[target]
        } else {
            bad.targets[target] <- new.targets[target]
        }
    }
    
    return(bad.targets)
}

MergeMisclassifications <- function(misclassifications, new.misclassifications) {
    for (target in names(new.misclassifications)) {
        if (target %in% names(misclassifications)) {
            misclassifications[[target]] <- append(misclassifications[[target]], new.misclassifications[[target]])
        } else {
            misclassifications[[target]] <- new.misclassifications[[target]]
        }
    }
    
    return(misclassifications)
}