UploadNewCompound <- function(feature.file.paths, similarity.file.paths, new.compound.features, new.compound.name){
    feature.file.paths <- list(sensitivity="Data/upload/features/sensitivity/aucs_all.RData", 
                               perturbation="Data/upload/features/perturbation/pert_features.RData",
                               structure="Data/upload/features/structure/structure_features.RData",
                               imaging=NULL)
    
    similarity.file.paths <- list(sensitivity="Data/upload/similarities/sensitivity/aucs_cor2.RData", 
                                  perturbation="Data/upload/similarities/perturbation/perturbation_similarities.RData",
                                  structure="Data/upload/similarities/structure/structure_similarities.RData",
                                  imaging=NULL)
    
    old.similarities <- lapply(similarity.file.paths, function(file.path) {
        if (!is.null(file.path)) {
            return(readRDS(file.path))
        } else {
            return(NULL)
        }
    })
    
    old.similarities[sapply(old.similarities, is.null)] <- NULL
    
    old.affinity.matrices <- ImputeSimilarities(old.similarities)
    integrated.old <- FuseAffinityMatrices(old.affinity.matrices)
    
    new.compound.name <- paste("user", new.compound.name, sep=".")
    
    res <- AddDrugToExistingSimilarities(feature.file.paths, similarity.file.paths, new.compound.features,
                                         new.compound.name)
    
    new.similarities <- res$similarities
    new.affinity.matrices <- ImputeSimilarities(new.similarities)
    integrated.new <- FuseAffinityMatrices(new.affinity.matrices)
    integrated.final <- FreezeNetwork(integrated.old, integrated.new)
    
    all.drugs <- colnames(integrated.final)
    data.bench <- DrugTargetsKNN(all.drugs, "temp", use.ctrpv2=TRUE,
                                 use.clue=TRUE, use.chembl=TRUE, use.dbank=TRUE, use.dtc=FALSE)
    nearest.neighbours <- GetKNearestNeighbours(new.compound.name, 7, integrated.final)
    targets <- GetTargetsForCompound(new.compound.name, nearest.neighbours, data.bench)
    
    position <- PositionDrugInNetwork(integrated)
    
    res <- list(integrated.final=integrated.final, targets=targets, nearest.neighbours=nearest.neighbours)
    
    return(res)
}

ImputeSimilarities <- function(correlation.matrices) {
    all.drugs <- c()
    
    for (i in names(correlation.matrices)) {
        all.drugs <- c(all.drugs, colnames(correlation.matrices[[i]]))
    }
    
    all.drugs <- sort(unique(all.drugs))
    
    correlation.matrices <- medianSimilarity(correlation.matrices)
    affinity.matrices <- CreateAffinityMatrices(correlation.matrices)
    augmented.matrices <- CreateAugmentedMatrixSkeletons(names(correlation.matrices), all.drugs)
    augmented.matrices <- ReplaceAugmentedExistingValues(augmented.matrices, affinity.matrices)
    affinity.matrices <- ReplaceAffinityMatrixValuesFast(augmented.matrices, correlation.matrices, 
                                                         all.drugs)
    affinity.matrices <- medianSimilarity(affinity.matrices)
    
    return(affinity.matrices)
}

FuseAffinityMatrices <- function(affinity.matrices, new.compound.name) {
    integrated <- SNFtool::SNF(affinity.matrices)
    rownames(integrated) <- rownames(affinity.matrices[[1]])
    colnames(integrated) <- colnames(affinity.matrices[[1]])
    
    return(integrated)
}

FreezeNetwork <- function(old.integrated, new.integrated) {
    frozen.network <- new.integrated
    frozen.network[rownames(old.integrated), colnames(old.integrated)] <- old.integrated[rownames(old.integrated), colnames(old.integrated)]
    
    return(frozen.network)
}

PositionDrugInNetwork <- function(integrated) {
    position <- list()
    
    return(position)
}

PredictTargetsForDrug <- function(integrated) {
    targets <- c() # named vector where values are probabilities
    
    return(targets)
}

GetKNearestNeighbours <- function(compound.name, k, integrated) {
    if (compound.name %in% rownames(integrated)) {
        all.neighbours <- integrated[compound.name, ]
    } else {
        # Throw exception here
    }
    
    if (k < 1) {
        # Throw excecption here
    }
    
    all.neighbours <- all.neighbours[!(names(all.neighbours) %in% compound.name)]
    k.nearest.neighbours <- head(sort(all.neighbours, decreasing = T), k)
    
    return(k.nearest.neighbours)
}

GetTargetsForCompound <- function(compound.name, nearest.neighbours, data.bench,
                                  scale.distance=FALSE, extra.scaling=100) {
    num.targets <- length(unique(data.bench$TARGET_NAME))
    drug.names <- sort(unique(data.bench$MOLECULE_NAME))
    target.names <- sort(unique(data.bench$TARGET_NAME))
    
    counts <- numeric(num.targets)
    names(counts) <- sort(unique(data.bench$TARGET_NAME))
    
    names.of.neighbours <- names(nearest.neighbours)
    
    max.x <- max(nearest.neighbours)
    fudge.factor <- - extra.scaling * ((log(3)) / ((nearest.neighbours[length(nearest.neighbours) - 1] - nearest.neighbours[length(nearest.neighbours)]))) 
    sigmoid <- function(z) {1 / (1 + exp(- fudge.factor * (z - max.x)))}
    
    for (j in 1:length(names.of.neighbours)) {
        neighbour.name <- names.of.neighbours[j]
        
        relevant.targets <- data.bench[data.bench$MOLECULE_NAME == neighbour.name, "TARGET_NAME"]
        
        if (scale.distance == TRUE) {
            counts[relevant.targets] <- counts[relevant.targets] + (sigmoid(nearest.neighbours[j]))    
        } else {
            counts[relevant.targets] <- counts[relevant.targets] + (nearest.neighbours[j])
        }
    }
    
    # Compute the softmax
    indices <- which(counts > 0)
    res <- numeric(length(counts))
    names(res) <- names(counts)
    res[indices] <- exp(counts[indices]) / sum(exp(counts[indices]))
    counts <- res    
    
    
    sorted <- sort(counts, index.return=T)
    index.of.predictions <- tail(sorted$ix[sorted$x > 0], 5)
    targets <- counts[index.of.predictions]
    
    return(targets)
}
