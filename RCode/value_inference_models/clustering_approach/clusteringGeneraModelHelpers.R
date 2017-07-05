GetNeighboursDistance <- function(features, drug, k=10) {
    neighbours <- c()
    
    distances <- rdist(features)
    rownames(distances) <- rownames(features)
    colnames(distances) <- rownames(features)
    relevant.distances <- distances[, drug]
    
    names(relevant.distances) <- rownames(distances)
    sorted.distances <- sort(relevant.distances)
    
    sorted.distances[1:k]
}

GetNeighbours <- function(precomputed.measure, drug, k=10) {
    neigbhours <- c()
    
    if (!(drug %in% colnames(precomputed.measure))) {
        return(c())
    }
    
    relevant.values <- precomputed.measure[, drug]
    names(relevant.values) <- rownames(precomputed.measure)
    
    #sorted.correlations <- sort(relevant.correlations, decreasing=TRUE)
    
    relevant.values
}

PredictValuesBasedOnNeighbours <- function(sens.data, pert.data, 
                                           imaging.data, drugs.to.predict) {
    precomputed.sens.measure <- NULL
    precomputed.pert.measure <- NULL
    
    if (!is.null(sens.data)) {
        #precomputed.sens.measure <- cor(t(sens.data), use="pairwise.complete.obs")
        precomputed.sens.measure <- rdist(sens.data)
        rownames(precomputed.sens.measure) <- rownames(sens.data)
        colnames(precomputed.sens.measure) <- rownames(sens.data)
    }
    
    if (!is.null(pert.data)) {
        #precomputed.pert.measure <- cor(t(pert.data), use="pairwise.complete.obs")
        precomputed.pert.measure <- rdist(pert.data)
        rownames(precomputed.pert.measure) <- rownames(pert.data)
        colnames(precomputed.pert.measure) <- rownames(pert.data)
    }
    
    for (drug in rownames(drugs.to.predict)) {
        sens.neighbours <- NULL
        pert.neighbours <- NULL
        
        if (!is.null(sens.data)) {
            sens.neighbours <- GetNeighbours(precomputed.sens.measure, drug, 200)    
        }
        
        if (!is.null(pert.data)) {
            pert.neighbours <- GetNeighbours(precomputed.pert.measure, drug, 200)
        }
        
        merged.neighbours <- Reduce(c, list(sens.neighbours, pert.neighbours))
        merged.neighbours <- merged.neighbours[names(merged.neighbours) %in% rownames(imaging.data)]
        merged.neighbours <- merged.neighbours[!duplicated(names(merged.neighbours))]
        
        merged.neighbours <- sort(merged.neighbours, decreasing=FALSE)[1:3]
        
        imaging.neighbours <- imaging.data[rownames(imaging.data) %in% names(merged.neighbours),]
        
        predicted.values <- colMeans(imaging.neighbours)
        drugs.to.predict[drug, ] <- predicted.values
    }
    
    drugs.to.predict
}