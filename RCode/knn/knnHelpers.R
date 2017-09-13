# This file contains functions that support the KNN target prediction.

GetKNNAccuracy <- function(k, integrated, data.bench, scale.distance=FALSE, extra.scaling=100, top=5) {
    # Creates a binary vector with an entry for each drug that indicates 
    # whether or not a single correct target was predicted to be in the top 
    # list of targets for a given drug.
    # 
    # Args:
    #   k: The number of neighbours to use in KNN.
    #   integrated: The SNF network to use in KNN.
    #   data.bench: A dataframe containing drugs and their targets.
    #   scale.distance: A logical value indicating whether or not to use a scaling 
    #                   function on network values.
    #   extra.scaling: A fudge factor to increase performance when using the scaling 
    #                  function on network values. Has no effect if scale.distance is FALSE.
    #   top: The number of top targets to consider for each drug. Any given drug could have 
    #        potentially tens of targets predicted for it, but this limits it to the top ones
    #        with highest confidence.
    #
    # Returns:
    #   A binary vector where each index corresponds to a certain drug and a value of 1 
    #   indicates that a true target was predicted for that drug, otherwise 0.
    drugs.with.targets <- sort(unique(data.bench$MOLECULE_NAME))
    integrated <- integrated[drugs.with.targets, drugs.with.targets]
    
    temp <- GetNearestNeighbours(k, integrated, drugs.with.targets)
    
    neighbours <- temp$neighbours
    weights <- temp$weights
    
    counts <- CountTargetsFromNeighbours(neighbours, weights, data.bench, 
                                         scale.distance, extra.scaling)
    top.predictions <- GetTopPredictions(counts, top)
    
    prediction.status <- EvaluatePredictionAccuracy(top.predictions, data.bench)
    
    return(prediction.status)
}

GetNearestNeighbours <- function(k, integrated, drugs) {
    # Determine k nearest neighbouurs for each drug in the provided integrated network.
    # 
    # Args:
    #   k: A numeric vector indicating the number of nearest neighbours to take into consideration.
    #   integrated: The SNF network to use in KNN.
    #   drugs: A character vector of drugs to subset to. This is usually the drugs for which there
    #          is drug target information available.
    #
    # Returns:
    #   A list containing two items. The first item in the list is a matrix that has the top k 
    #   neighbours for each drug in the argument drugs provided. The second element in the list 
    #   is also a matrix that indicates the similarity of a given drug with its various k nearest
    #   neighbours.
    drugs.with.targets <- rownames(integrated)
    
    if (k == ncol(integrated)) {
        k <- k - 1
    }
    
    # Neighbours is a a matirx of the shape NUM_DRUGS X K. 
    # For example, a row looks like: Drug X     34 21 ... 10 4.
    neighbours <- matrix(0, nrow=length(drugs), ncol=k)
    rownames(neighbours) <- drugs
    
    weights <- matrix(0, nrow=length(drugs), ncol=k)
    rownames(weights) <- drugs
    
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
        
        res.indices <- res.indices[!(names(res.indices) %in% drug.name)]
        res.weights <- res.weights[!(names(res.weights) %in% drug.name)]

        res.indices <- res.indices[1:(min(k, length(res.indices)))]
        res.weights <- res.weights[1:(min(k, length(res.weights)))]

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
        rownames(neighbours) <- drugs        
    }

    return(list(neighbours=neighbours, weights=weights))
}

CountTargetsFromNeighbours <- function(neighbours, weights, data.bench,
                                       scale.distance=FALSE, extra.scaling=100) {
    
    # Creates a matrix of accumulated similarties. More specifically, each row of the matrix  
    # corresponds to a drug, and each column corresponds to a possible target in our drug target
    # list. Entries in the matrix accumulate values as follows. For each neighbour of a given drug X,
    # determine that neighbour's true targets from the drug target bencmark. Then, use those targets
    # to index the matrix of accumulated similarities and add the similarity between drug X and its neighbour
    # to all cells in row X. This happens for all drugs.
    #
    # Args:
    #   neighbours: An N x K character matrix (N is the total number of drugs, K is the number of nearest neighbours)
    #               where rownames are drugs, and each row contains the k nearest 
    #               neighbours of the drug specified in the rowname.
    #   weights: An N x K numeric matrix (N is the total number of drugs, K is the number of nearest neighbours)
    #            where rownames are drugs, and each row contains the simliarities of the k nearest neighbours
    #            with the drug specified in the rowname.
    #   data.bench: A dataframe containing containing two columns: one column is molecule name, and
    #               the other column is target name. This is the drug target ground truth.
    #   scale.distance: A logical value indicating whether or not to use a scaling 
    #                   function on network values.
    #   extra.scaling: A fudge factor to increase performance when using the scaling 
    #                  function on network values. Has no effect if scale.distance is FALSE.
    #
    # Returns: 
    #   A matrix of counts (usually fractional counts) where rows correspond to drugs and columns 
    #   correspond to targets. The value of a given cell (i, j) in the matrix indicates the likelihood (abuse of 
    #   probability notation) that a given drug i has target j.
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
    # Given a matrix of drug target counts, return a list containing the top targets
    # for each drug according to count.
    # 
    # Args:
    #   counts: A numeric matrix where each row corresponds to a drug and each column corresponds to a target.
    #           The value of a given cell (i, j) in the matrix indicates the likelihood (abuse of 
    #           probability notation) that a given drug i has target j.
    #   top: The number of top targets to consider for each drug. Any given drug could have 
    #        potentially tens of targets predicted for it, but this limits it to the top ones
    #        with highest confidence.
    # 
    # Returns:
    #   A list of the top most likely targets for each drug.
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
    # Creates a binary vector indicating where each entry corresponds to a drug and takes
    # on the value of 1 if a target was correctly predicted for that drug, and 0 otherwise.
    #
    # Args:
    #   top.predictions: A list of drug target predictions.
    #   data.bench: A dataframe containing containing two columns: one column is molecule name, and
    #               the other column is target name. This is the drug target ground truth.
    #
    # Returns:
    #   A binary vector indicating whether or not a correct drug target prediction was made for each drug.
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

GetKNNAccuracySingle <- function(k, integrated, data.bench, drug, scale.distance=FALSE, extra.scaling=100, top=5) {
    # Makes use of the KNN drug target prediction method to evaluate single-drug performance.
    #
    # Args:
    #   k: The number of neighbours to use in KNN.
    #   integrated: The SNF network to use in KNN.
    #   data.bench: A dataframe containing drugs and their targets.
    #   scale.distance: A logical value indicating whether or not to use a scaling 
    #                   function on network values.
    #   extra.scaling: A fudge factor to increase performance when using the scaling 
    #                  function on network values. Has no effect if scale.distance is FALSE.
    #   top: The number of top targets to consider for each drug. Any given drug could have 
    #        potentially tens of targets predicted for it, but this limits it to the top ones
    #        with highest confidence.
    #
    # Returns:
    #   A binary value taking on 1 to indicate target prediction success, 0 otherwise.
    
    drugs.with.targets <- sort(unique(data.bench$MOLECULE_NAME))
    integrated <- integrated[drugs.with.targets, drugs.with.targets]
    
    
    temp <- GetNearestNeighbours(k, integrated, drug)
    
    neighbours <- temp$neighbours
    weights <- temp$weights
    
    counts <- CountTargetsFromNeighbours(neighbours, weights, data.bench, 
                                         scale.distance, extra.scaling)
    top.predictions <- GetTopPredictions(counts, top)
    
    prediction.status <- EvaluatePredictionAccuracy(top.predictions, data.bench)
    
    return(prediction.status)
}
