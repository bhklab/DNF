GetKNNAccuracy <- function(k, integrated, data.bench, scale.distance=FALSE, extra.scaling=100, top=5) {
    drugs.with.targets <- sort(unique(data.bench$MOLECULE_NAME))
    
    temp <- GetNearestNeighbours(k, integrated, drugs.with.targets)
    
    neighbours <- temp$neighbours
    weights <- temp$weights
    
    counts <- CountTargetsFromNeighbours(neighbours, weights, data.bench, 
                                         scale.distance, extra.scaling)
    top.predictions <- GetTopPredictions(counts, top)
    
    prediction.status <- EvaluatePredictionAccuracy(top.predictions, data.bench)
    
    return(prediction.status)
}

GetNearestNeighbours <- function(k, integrated, drugs.with.targets) {
    integrated <- integrated[drugs.with.targets, drugs.with.targets]
    
    # Neighbours is a a matirx of the shape NUM_DRUGS X K. 
    # For example, a row looks like: Drug X     34 21 ... 10 4.
    neighbours <- matrix(0, nrow=length(drugs.with.targets), ncol=k)
    rownames(neighbours) <- drugs.with.targets
    
    weights <- matrix(0, nrow=length(drugs.with.targets), ncol=k)
    rownames(weights) <- drugs.with.targets
    
    # For each drug, obtain the indices of the k closest neighbours
    # according to the integrated similarity measure. For some reason,
    # using existing KNN packages gave different results (I guess they
    # thought that rows corresponds to points and columns correspond to 
    # dimensions).
    for (i in 1:nrow(neighbours)) {
        drug.name <- rownames(neighbours)[i]
        row <- integrated[drug.name, ]
        
        res <- sapply(sort(row, index.return=TRUE, decreasing=T), '[')
        res.indices <- res[, "ix"]
        res.weights <- res[, 'x']

        res.indices <- res.indices[2:(k+1)]
        res.weights <- res.weights[2:(k+1)]
        
        res.indices <- rev(res.indices)
        res.weights <- rev(res.weights)

        neighbours[i, ] <- res.indices
        weights[i, ] <- res.weights
    }
    
    # Convert the indices of nearest neigbhours into the actual drug
    # names of those neighbours.
    neighbours <- apply(neighbours, 1, function(x) {drugs.with.targets[x]})
    if (k == 1) {
        neighbours <- as.matrix(neighbours)    
    } else {
        neighbours <- t(neighbours)
        rownames(neighbours) <- drugs.with.targets        
    }

    return(list(neighbours=neighbours, weights=weights))
}

CountTargetsFromNeighbours <- function(neighbours, weights, data.bench,
                                       scale.distance=FALSE, extra.scaling=100) {
    num.targets <- length(unique(data.bench$TARGET_NAME))
    counts <- matrix(data=0, nrow=nrow(neighbours), ncol=num.targets)
    rownames(counts) <- rownames(neighbours)
    colnames(counts) <- sort(unique(data.bench$TARGET_NAME))
    
    # Now if we want to weight the targets based on the similarities between drugs,
    # we have to determine some kind of threshold. For example, if the similarity
    # between drug X and drug Y is 0.003, and the similarity between drug X and drug
    # Z is 0.001, we should definitely trust the targets coming from drug Y more than 
    # those coming from drug Z.
    
    # Iterate over each drug in the neighbours matrix. For each neighbour
    # of a given drug, obtain that neighbour's drug targets and use those 
    # targets to increment the counts matrix in the apporpriate location.
    for (i in 1:nrow(neighbours)) {
        drug <- rownames(neighbours)[i]
        names.of.neighbours <- neighbours[i, ]
        
        max.x <- max(weights[i, ])
        fudge.factor <- - extra.scaling * ((log(3)) / ((weights[i, ncol(weights) - 1] - weights[i, ncol(weights)]))) 
        sigmoid <- function(z) {1 / (1 + exp(- fudge.factor * (z - max.x)))}
        
        for (j in 1:length(names.of.neighbours)) {
            neighbour.name <- names.of.neighbours[j]
            
            relevant.targets <- data.bench[data.bench$MOLECULE_NAME == neighbour.name, "TARGET_NAME"]
            
            if (scale.distance == TRUE & ncol(counts) > 1) {
                counts[drug, relevant.targets] <- counts[drug, relevant.targets] + (sigmoid(weights[i, j]))    
            } else {
                counts[drug, relevant.targets] <- counts[drug, relevant.targets] + (weights[i, j])
            }
        }
    }
    
    return(counts)
}

GetTopPredictions <- function(counts, top=5) {
    top.predictions <- sapply(1:nrow(counts), function(i, target.names, drugs.with.targets) {
        x <- counts[i, ]
        sorted <- sort(x, index.return=T, decreasing=T)
        index.of.predictions <- head(sorted$ix[sorted$x > 0], top)
        res <- list()
        res[[drugs.with.targets[i]]] <- x[index.of.predictions]
        names(res[[drugs.with.targets[i]]]) <- target.names[index.of.predictions]
        res
    }, target.names=colnames(counts), drugs.with.targets=rownames(counts))
    
    return(top.predictions)
}

EvaluatePredictionAccuracy <- function(top.predictions, data.bench) {
    prediction.status <- integer(length(top.predictions))
    names(prediction.status) <- names(top.predictions)
    
    for (i in 1:length(top.predictions)) {
        drug.name <- names(top.predictions)[i]
        target.names <- names(top.predictions[[i]])
        
        possible.true.targets <- data.bench[data.bench$MOLECULE_NAME == drug.name, "TARGET_NAME"]
        
        if (length(intersect(target.names, possible.true.targets)) > 0) {
            prediction.status[drug.name] <- 1
        }
    }
    
    return(prediction.status)
}

GetKNNAccuracySingleDrug <- function(k, integrated, data.bench, drug.name, scale.distance=FALSE, extra.scaling=100) {
    num.targets <- length(unique(data.bench$TARGET_NAME))
    drugs.with.targets <- sort(unique(data.bench$MOLECULE_NAME))
    target.names <- sort(unique(data.bench$TARGET_NAME))
    
    counts <- matrix(data=0, nrow=length(drugs.with.targets), ncol=num.targets)
    rownames(counts) <- drugs.with.targets
    colnames(counts) <- sort(unique(data.bench$TARGET_NAME))
    
    # Neighbours is a a matirx of the shape NUM_DRUGS X K. 
    # For example, a row looks like: Drug X     34 21 ... 10 4.
    neighbours <- matrix(0, nrow=1, ncol=k)
    rownames(neighbours) <- drug.name
    
    weights <- matrix(0, nrow=1, ncol=k)
    rownames(weights) <- drug.name
    
    # For each drug, obtain the indices of the k closest neighbours
    # according to the integrated similarity measure. For some reason,
    # using existing KNN packages gave different results (I guess they
    # thought that rows corresponds to points and columns correspond to 
    # dimensions).
    row <- integrated[drug.name, ]
    res <- sapply(sort(row, index.return=TRUE), '[')
    res.indices <- res[, "ix"]
    res.weights <- res[, 'x']
    
    res.indices <- res.indices[(length(res.indices) - (k)):(length(res.indices) - 1)]
    res.weights <- res.weights[(length(res.weights) - (k)):(length(res.weights) - 1)]
    neighbours[drug.name, ] <- res.indices
    weights[drug.name, ] <- res.weights
    
    # Convert the indices of nearest neigbhours into the actual drug
    # names of those neighbours.
    neighbours <- apply(neighbours, 1, function(x) {drugs.with.targets[x]})
    neighbours <- t(neighbours)
    rownames(neighbours) <- drug.name
    
    # Now if we want to weight the targets based on the similarities between drugs,
    # we have to determine some kind of threshold. For example, if the similarity
    # between drug X and drug Y is 0.003, and the similarity between drug X and drug
    # Z is 0.001, we should definitely trust the targets coming from drug Y more than 
    # those coming from drug Z.
    
    # Iterate over each drug in the neighbours matrix. For each neighbour
    # of a given drug, obtain that neighbour's drug targets and use those 
    # targets to increment the counts matrix in the apporpriate location.
    drug <- rownames(neighbours)[1]
    names.of.neighbours <- neighbours[1, ]
    
    max.x <- max(weights[1, ])
    fudge.factor <- - extra.scaling * ((log(3)) / ((weights[1, ncol(weights) - 1] - weights[1, ncol(weights)]))) 
    sigmoid <- function(z) {1 / (1 + exp(- fudge.factor * (z - max.x)))}
    
    for (j in 1:length(names.of.neighbours)) {
        neighbour.name <- names.of.neighbours[j]
        
        relevant.targets <- data.bench[data.bench$MOLECULE_NAME == neighbour.name, "TARGET_NAME"]
        
        if (scale.distance == TRUE) {
            counts[drug, relevant.targets] <- counts[drug, relevant.targets] + (sigmoid(weights[1, j]))    
        } else {
            counts[drug, relevant.targets] <- counts[drug, relevant.targets] + (weights[1, j])
        }
    }
    
    counts <- apply(counts, 1, function(x) {
        indices <- which(x > 0)
        res <- numeric(length(x))
        res[indices] <- exp(x[indices]) / sum(exp(x[indices]))
        res
    })
    
    counts <- t(counts)
    colnames(counts) <- sort(unique(data.bench$TARGET_NAME))
    
    top.predictions <- sapply(1:nrow(counts), function(i, target.names, drugs.with.targets) {
        x <- counts[i, ]
        sorted <- sort(x, index.return=T)
        index.of.predictions <- tail(sorted$ix[sorted$x > 0], 5)
        res <- list()
        res[[drugs.with.targets[i]]] <- x[index.of.predictions]
        names(res[[drugs.with.targets[i]]]) <- target.names[index.of.predictions]
        res
    }, target.names=colnames(counts), drugs.with.targets=rownames(counts))
    
    for (i in 1:length(top.predictions)) {
        drug.name <- names(top.predictions)[i]
        target.names <- names(top.predictions[[i]])
        
        possible.true.targets <- data.bench[data.bench$MOLECULE_NAME == drug.name, "TARGET_NAME"]
        
        if (length(intersect(target.names, possible.true.targets)) > 0) {
            num.correct.predictions[drug.name] <- 1
        }
    }
    
    accuracy <- sum(num.correct.predictions) / length(num.correct.predictions)
    
    accuracy
}