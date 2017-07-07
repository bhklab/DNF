# This file contains helper functions that can add new features to existing similarity matrices
# I think we should separate functionality as much as we possibly can.
# There should be a function that just loads in features from files
# There should be another function that 

LoadSensitivityFeatures <- function(file.path) {
    sensitivity.features <- readRDS(file.path)
    
    sensitivity.features
}

LoadSensitivitySimilarities <- function(file.path) {
    sensitivity.similarities <- readRDS(file.path)
    
    sensitivity.similarities
}

AddNewSensitivityFeatures <- function(existing.features, existing.similarities, new.compound.name, new.compound.features) {
    concatenated.drug.name <- paste("new.dataset", new.compound.name, sep=":::")
    
    if (new.compound.name %in% colnames(existing.features)) {
        return(existing.features)
    }
    
    new.col <- rep(NA, nrow(existing.features))
    existing.features <- cbind(existing.features, new.col)
    
    colnames(existing.features)[ncol(existing.features)] <- concatenated.drug.name
    
    
    # Now we can set the values of the new row according to new.compound.features
    common.feature.names <- intersect(rownames(existing.features), names(new.compound.features))
    existing.features[common.feature.names, concatenated.drug.name] <- new.compound.features[common.feature.names] 
    
    new.similarities <- ComputeSensitivityCorrelations(existing.features, 
                                                       new.compound.features,
                                                       common.feature.names)
    
    existing.similarities <- rbind(existing.similarities, numeric(ncol(existing.similarities)))
    existing.similarities <- cbind(existing.similarities, numeric(nrow(existing.similarities)))
    
    colnames(existing.similarities)[ncol(existing.similarities)] <- new.compound.name
    rownames(existing.similarities)[nrow(existing.similarities)] <- new.compound.name
    
    existing.similarities[new.compound.name, ] <- new.similarities[rownames(existing.similarities)]
    existing.similarities[, new.compound.name] <- new.similarities[colnames(existing.similarities)]
    existing.similarities <- existing.similarities[order(rownames(existing.similarities)), 
                                                   order(colnames(existing.similarities))]
    
    list(updated.features=existing.features, updated.similarities=existing.similarities)
}

ComputeSensitivityCorrelations <- function(existing.features, new.compound.features, common.feature.names) {
    new.similarities <- apply(existing.features, 2, cor, y = new.compound.features[common.feature.names],
                              use="pairwise.complete.obs")
    names(new.similarities) <- colnames(existing.features)
    
    num.samples.used <- apply(existing.features, 2, function(x, y) {
        x.len <- length(x[!is.na(x)])
        y.len <- length(y[!is.na(y)])
        
        min(x.len, y.len)
    }, y=new.compound.features)
    
    names(num.samples.used) <- colnames(existing.features)
    
    fisher.transformed <- fisherz(new.similarities)
    
    
    standard.errors <- 1 / sqrt(num.samples.used - 3)
    names(standard.errors) <- names(num.samples.used)
    
    unique.drug.names <- unique(sapply(strsplit(names(new.similarities), ":::"), function(x) {return(x[2])}))
    final.correlations <- numeric(length(unique.drug.names))
    names(final.correlations) <- unique.drug.names
    
    for (i in 1:length(unique.drug.names)) {
        drug.name <- unique.drug.names[i]
        
        duplicated.drugs <- colnames(existing.features)[endsWith(colnames(existing.features),
                                                                 paste(":", drug.name, sep=""))]
        
        combined <- combine.est(fisher.transformed[duplicated.drugs], num.samples.used[duplicated.drugs])
        final.correlations[drug.name] <- fisherz(combined$estimate, inv=TRUE)
    }
    
    final.correlations
}

LoadStructureFeatures <- function(file.path) {
    structure.features <- readRDS(file.path)
    
    structure.features
}

LoadStructureSimilarities <- function(file.path) {
    structure.similarities <- readRDS(file.path)
    
    structure.similarities
}

AddNewStructureFeatures <- function(existing.features, existing.similarities, new.compound.name, new.compound.smiles) {
    if (new.compound.name %in% names(existing.features)) {
        res <- list(updated.features=existing.features, updated.similarities=existing.similarities)
        return(res)
    }
    
    new.fingerprint <- rcdk::parse.smiles(new.compound.smiles)
    new.fingerprint <- rcdk::get.fingerprint(new.fingerprint)
    existing.features[[new.compound.name]] <- new.fingerprint
    
    new.similarities <- numeric(ncol(existing.similarities))
    new.similarities <- apply(existing.features, 1, cor, y = new.compound.features[common.feature.names],
                              use="pairwise.complete.obs")
    new.similarities <- c(new.similarities, 1)
    
    existing.similarities <- rbind(existing.similarities, numeric(ncol(existing.similarities)))
    existing.similarities <- cbind(existing.similarities, numeric(nrow(existing.similarities)))
    
    colnames(existing.similarities)[ncol(existing.similarities)] <- new.compound.name
    rownames(existing.similarities)[nrow(existing.similarities)] <- new.compound.name
    
    existing.similarities[new.compound.name, ] <- new.similarities
    existing.similarities[, new.compound.name] <- new.similarities
    existing.similarities <- existing.similarities[order(rownames(existing.similarities)), 
                                                   order(colnames(existing.similarities))]
    
    list(updated.features=existing.features, updated.similarities=existing.similarities)
}

LoadPerturbationFeatures <- function(file.path) {
    perturbation.features <- readRDS(file.path)
    
    perturbation.features
}

LoadPerturbationSimilarities <- function(file.path) {
    perturbation.similarities <- readRDS(file.path)
    
    perturbation.similarities
}

AddNewPerturbationFeatures <- function(existing.features, existing.similarities, new.compound.name, new.compound.features) {
    # new.compound.features should be a named vector of some sort
    # existing.features should be a matrix of the form NUM_DRUGS X NUM_FEATURES
    
    # Let's first check that the new compound is in fact new
    if (new.compound.name %in% colnames(existing.features)) {
        res <- list(updated.features=existing.features, updated.similarities=existing.similarities)
        return(res)
    }
    
    # Now let's add a new row to existing.features
    new.col <- rep(NA, nrow(existing.features))
    existing.features <- cbind(existing.features, new.col)
    colnames(existing.features)[ncol(existing.features)] <- new.compound.name
    
    # Now we can set the values of the new row according to new.compound.features
    common.feature.names <- intersect(rownames(existing.features), names(new.compound.features))
    existing.features[common.feature.names, new.compound.name] <- new.compound.features[common.feature.names] 

    new.similarities <- numeric(nrow(existing.similarities))
    new.similarities <- apply(existing.features, 2, cor, y = new.compound.features[common.feature.names],
                              use="pairwise.complete.obs")
    # new.similarities <- c(new.similarities, 1)
    
    existing.similarities <- rbind(existing.similarities, numeric(ncol(existing.similarities)))
    existing.similarities <- cbind(existing.similarities, numeric(nrow(existing.similarities)))
    
    colnames(existing.similarities)[ncol(existing.similarities)] <- new.compound.name
    rownames(existing.similarities)[nrow(existing.similarities)] <- new.compound.name
    
    existing.similarities[new.compound.name, ] <- new.similarities[rownames(existing.similarities)]
    existing.similarities[, new.compound.name] <- new.similarities[colnames(existing.similarities)]
    existing.similarities <- existing.similarities[order(rownames(existing.similarities)), 
                                                   order(colnames(existing.similarities))]
    
    list(updated.features=existing.features, updated.similarities=existing.similarities)
}

LoadImagingFeatures <- function(file.path) {
    imaging.features <- readRDS(file.path)
    
    imaging.features
}

LoadImagingSimilarities <- function(file.path) {
    imaging.similarities <- readRDS(file.path)
    
    imaging.similarities
}

AddNewImagingFeatures <- function(existing.features, existing.similarities, new.compound.name, new.compound.features) {
    if (new.compound.name %in% rownames(existing.features)) {
        return(existing.features)
    }
    
    new.row <- numeric(ncol(existing.features))
    existing.features <- rbind(existing.features, new.row)
    rownames(existing.features)[nrow(existing.features)] <- new.compound.name
    
    common.feature.names <- intersect(colnames(existing.features), names(new.compound.features))
    existing.features[new.compound.name, common.feature.names] <- new.compound.features[common.feature.names]
    
    existing.features
}

UploadFeaturesForNewCompound <-  function(new.features, new.compound.name) {
    # Pretend like there is some global list of file names
    # features is a list where the name of each list element is the layer type
    # and the values are the features for the new compound corresponding to 
    # the said layer type

    structure.features <- NULL
    perturbation.features <- NULL
    sensitivity.features <- NULL
    imaging.features <- NULL
    
    #structure.features <- LoadStructureFeatures(feature.file.paths[['structure']])
    perturbation.features <- LoadPerturbationFeatures(feature.file.paths[['perturbation']])
    sensitivity.features <- LoadSensitivityFeatures(feature.file.paths[['sensitivity']])
    #imaging.features <- LoadImagingFeatures(feature.file.paths[['imaging']])

    existing.features <- list(sensitivity=sensitivity.features, perturbation=perturbation.features,
                              structure=structure.features, imaging=imaging.features)
    
    structure.similarities <- NULL
    perturbation.similarities <- NULL
    sensitivity.similarities <- NULL
    imaging.similarities <- NULL
        
    #structure.similarities <- LoadStructureSimilarities(similarity.file.paths[['structure']])
    perturbation.similarities <- LoadPerturbationSimilarities(similarity.file.paths[['perturbation']])
    sensitivity.similarities <- LoadSensitivitySimilarities(similarity.file.paths[['sensitivity']])
    # imaging.similarities <- LoadImagingSimilarities(similarity.file.paths[['imaging']])
    
    existing.similarities <- list(sensitivity=sensitivity.similarities, perturbation=perturbation.similarities,
                                  structure=structure.similarities, imaging=imaging.similarities)
    
    # We iterate over the features for the new drug and add the features to the 
    # apporpriate layers. Then, for the layers that we didn't add features to,
    # we infer the features for those layers.
    
    for (layer in names(new.features)) {
        if (is.null(new.features[[layer]])) {
            next()
        }
        
        new.compound.features <- new.features[[layer]]
        
        if (layer == 'sensitivity') {
            res <- AddNewSensitivityFeatures(existing.features[[layer]], existing.similarities[[layer]],
                                             new.compound.name, new.compound.features)
            
            existing.features[[layer]] <- res$updated.features
            existing.similarities[[layer]] <- res$updated.similarities
            
        } else if (layer == 'perturbation') {
            res <- AddNewPerturbationFeatures(existing.features[[layer]], existing.similarities[[layer]],
                                              new.compound.name, new.compound.features)
            
            existing.features[[layer]] <- res$updated.features
            existing.similarities[[layer]] <- res$updated.similarities
        } else if (layer == 'structure') {
            res <- AddNewStructureFeatures(existing.features[[layer]], existing.similarities[[layer]],
                                           new.compound.name, new.compound.features)
            
            existing.features[[layer]] <- res$updated.features
            existing.similarities[[layer]] <- res$updated.similarities
        } else if (layer == 'imaging') {
            res <- AddNewImagingFeatures(existing.features[[layer]], existing.similarities[[layer]],
                                         new.compound.name, new.compound.features)
            
            existing.features[[layer]] <- res$updated.features
            existing.similarities[[layer]] <- res$updated.similarities
        }
    }
    
    # Now for the layers that the user didn't upload data for, we can use our predictive models
    # to infer the values for these layers. We should never have to infer values for the structure
    # layer however. Structure information should always be available for a compound.
    # For now, since we don't have all the necessary models built yet, we can just impute using
    # missForest. However, this approach probably will give worse results than a predictive model.
    # Furthermore, imputation can't be used for more than one drug. This is because if we have 2 or
    # more NA rows, those drugs will probably end up getting the same imputed values, so they will have
    # exactly the same features, and end getting a drug-drug similarity score of 1! Maybe the clustering
    # approach would be best, though it remains to be determined.
    
    # c1 <- makeCluster(8)
    # registerDoParallel(c1)
    # for (layer in names(new.features)) {
    #     if (!is.null(new.features[[layer]])) {
    #         next()
    #     }
    # 
    #     # Like we said above, for now we just impute data until we have some proper predictive models.
    # 
    #     new.row <- numeric(ncol(existing.features[[layer]]))
    #     existing.features[[layer]] <- rbind(existing.features[[layer]], new.row)
    #     rownames(existing.features[[layer]])[nrow(existing.features[[layer]])] <- new.compound.name
    #     existing.features[[layer]] <- missForest::missForest(existing.features[[layer]], maxiter=10,
    #                                                          parallelize = 'variables')
    # }
}