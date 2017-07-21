# This file contains helper functions that can add new features to existing similarity matrices
# I think we should separate functionality as much as we possibly can.
# There should be a function that just loads in features from files
# There should be another function that 
library(survcomp)

LoadSensitivityFeatures <- function(file.path) {
    # Load sensitivity features from the specified file path.
    #
    # Args:
    #   file.path A character(1) containing the file path.
    #
    # Returns:
    #   Sensitivity feaures from the specified file path.
    sensitivity.features <- readRDS(file.path)
    
    sensitivity.features
}

LoadSensitivitySimilarities <- function(file.path) {
    # Load sensitivity similarities from the specified file path.
    #
    # Args:
    #   file.path A character(1) containing the file path.
    #
    # Returns:
    #   Sensitivity similarities from the specified file path.
    sensitivity.similarities <- readRDS(file.path)
    
    sensitivity.similarities
}

AddNewSensitivityFeatures <- function(existing.features, existing.similarities, new.compound.name, new.compound.features) {
    # Add sensitivity features for a new compound to a matrix of already existing sensitivity features. Also,
    # augment the existing sensitivity similarity matrix to include the similarities with the new compound.
    #
    # Args:
    #   existing.features: A matrix of existing sensitivity features where rownames correspond to cell
    #   lines and colnames contain drugnames preprended with the sensitivity dataset name.
    #   existing.similarities: A correlation matrix of existing similarities. Rownames and colnames are both
    #   are both just drug names.
    #   new.compound.name: A character(1) for the name of the new compound.
    #   new.compound.features: A named numeric vector containing the sensitivity features for the new 
    #   compound.
    #
    # Returns:
    #   A list with two elements. The first element is an updated feature matrix, and the 
    #   second element is an updated similarity matrix.
    concatenated.drug.name <- paste("new.dataset", new.compound.name, sep=":::")
    
    # If the new.compound.name isn't actually new, then return existing features and similarities.
    if (new.compound.name %in% colnames(existing.similarities)) {
        return(list(updated.features=existing.features, updated.similarities=existing.similarities))
    }
    
    # If there are no features in common between the new compound and existing features,
    # then return the existing features and similarities.
    if (length(intersect(names(new.compound.features), rownames(existing.features))) == 0) {
        return(list(updated.features=existing.features, updated.similarities=existing.similarities))
    }
    
    # Add a new column to the feature matrix to hold the values of the features for 
    # the new compound.
    new.col <- rep(NA, nrow(existing.features))
    existing.features <- cbind(existing.features, new.col)
    colnames(existing.features)[ncol(existing.features)] <- concatenated.drug.name
    
    
    # Now we can set the values of the new column according to new.compound.features
    common.feature.names <- intersect(rownames(existing.features), names(new.compound.features))
    existing.features[common.feature.names, concatenated.drug.name] <- new.compound.features[common.feature.names] 
    
    # Compute correlation of the new compound with existing compounds.
    new.similarities <- ComputeSensitivityCorrelations(existing.features, 
                                                       new.compound.features,
                                                       common.feature.names)
    
    # Augment the similarity matrix to include correlations with the new compound.
    existing.similarities <- rbind(existing.similarities, numeric(ncol(existing.similarities)))
    existing.similarities <- cbind(existing.similarities, numeric(nrow(existing.similarities)))
    
    colnames(existing.similarities)[ncol(existing.similarities)] <- new.compound.name
    rownames(existing.similarities)[nrow(existing.similarities)] <- new.compound.name
    
    existing.similarities[new.compound.name, ] <- new.similarities[colnames(existing.similarities)]
    existing.similarities[, new.compound.name] <- new.similarities[rownames(existing.similarities)]
    existing.similarities <- existing.similarities[order(rownames(existing.similarities)), 
                                                   order(colnames(existing.similarities))]
    
    list(updated.features=existing.features, updated.similarities=existing.similarities)
}

ComputeSensitivityCorrelations <- function(existing.features, new.compound.features, common.feature.names) {
    placeholder <- rep(NA, nrow(existing.features))
    names(placeholder) <- rownames(existing.features)
    
    placeholder[common.feature.names] <- new.compound.features[common.feature.names]
    
    new.similarities <- apply(existing.features, 2, cor, y = placeholder,
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
    # Load structure features from the specified file path.
    #
    # Args:
    #   file.path A character(1) containing the file path.
    #
    # Returns:
    #   Structure feaures from the specified file path.
    structure.features <- readRDS(file.path)
    
    structure.features
}

LoadStructureSimilarities <- function(file.path) {
    # Load structure similarities from the specified file path.
    #
    # Args:
    #   file.path A character(1) containing the file path.
    #
    # Returns:
    #   Structure similarities from the specified file path.
    structure.similarities <- readRDS(file.path)
    
    structure.similarities
}

AddNewStructureFeatures <- function(existing.features, existing.similarities, new.compound.name, new.compound.smiles) {
    # Add structure features for a new compound to a list of already existing structure features. Also,
    # augment the existing structure similarity matrix to include similarities with the new compound.
    #
    # Args:
    #   exiting.features: A list of chemical fingerprints obtained with the fingerprint package.
    #   existing.similarities: A similarity matrix (Tanimoto) for the existing features. Rownames
    #   and colnames are both just drug names.
    #   new.compound.name: A character(1) for the name of the new compound.
    #   new.compound.smiles: A character(1) containing the SMILES representation of the new compound
    #
    # Returns:
    #   A list with two elements. The first element is an updated list of fingerprints, and the 
    #   second element is an updated similarity matrix.
    
    # If the new compound isn't actually new, then return the existing features and similarities.
    if (new.compound.name %in% names(existing.features)) {
        res <- list(updated.features=existing.features, updated.similarities=existing.similarities)
        return(res)
    }
    
    # Create a fingerprint for the new compound based on the provided SMILES string.
    new.fingerprint <- rcdk::parse.smiles(new.compound.smiles)
    
    # If the parsed SMILES is NA, then the SMILES that the user provided is invalid. For now,
    # return the existing features and similarities.
    if (is.na(new.fingerprint)) {
        res <- list(updated.features=existing.features, updated.similarities=existing.similarities)
        return(res)
    }
    
    new.fingerprint <- rcdk::get.fingerprint(new.fingerprint[[1]], type="extended")
    
    # Add the new fingerprint to the list of existing fingerprints
    existing.features[[new.compound.name]] <- new.fingerprint
    
    # Calculate the similarities of the new fingerprint with the existing fingerprints.
    new.similarities <- unlist(lapply(existing.features, fingerprint::distance, fp2=new.fingerprint))
    
    # Augment the similarity matrix to include the similarities with the new compound.
    existing.similarities <- rbind(existing.similarities, numeric(ncol(existing.similarities)))
    existing.similarities <- cbind(existing.similarities, numeric(nrow(existing.similarities)))
    
    colnames(existing.similarities)[ncol(existing.similarities)] <- new.compound.name
    rownames(existing.similarities)[nrow(existing.similarities)] <- new.compound.name
    
    existing.similarities[new.compound.name, ] <- new.similarities[colnames(existing.similarities)]
    existing.similarities[, new.compound.name] <- new.similarities[rownames(existing.similarities)]
    existing.similarities <- existing.similarities[order(rownames(existing.similarities)), 
                                                   order(colnames(existing.similarities))]
    
    list(updated.features=existing.features, updated.similarities=existing.similarities)
}

LoadPerturbationFeatures <- function(file.path) {
    # Load perturbation features from the specified file path.
    #
    # Args:
    #   file.path A character(1) containing the file path.
    #
    # Returns:
    #   Perturbation feaures from the specified file path.
    perturbation.features <- readRDS(file.path)
    
    perturbation.features
}

LoadPerturbationSimilarities <- function(file.path) {
    # Load perturbation similarities from the specified file path.
    #
    # Args:
    #   file.path A character(1) containing the file path.
    #
    # Returns:
    #   Perturbation similarities from the specified file path.
    perturbation.similarities <- readRDS(file.path)
    
    perturbation.similarities
}

AddNewPerturbationFeatures <- function(existing.features, existing.similarities, new.compound.name, new.compound.features) {
    # Add pertrubation features for a new compound to a matrix of already existing perturbation features. Also,
    # augment the existing perturbation similarity matrix to include the similarities with the new compound.
    #
    # Args:
    #   existing.features: A matrix of existing perturbation features where rownames correspond to genes
    #   and colnames contain drug names.
    #   existing.similarities: A correlation matrix of existing perturbation similarities. Rownames and 
    #   colnames are both are both just drug names.
    #   new.compound.name: A character(1) for the name of the new compound.
    #   new.compound.features: A named numeric vector containing the perturbation features for the new 
    #   compound.
    #
    # Returns:
    #   A list with two elements. The first element is an updated feature matrix, and the 
    #   second element is an updated similarity matrix.
    
    # If the new.compound.name isn't actually new, then return existing features and similarities.
    if (new.compound.name %in% colnames(existing.features)) {
        return(list(updated.features=existing.features, updated.similarities=existing.similarities))
    }
    
    # If there are no features in common between the new compound and existing features,
    # then return the existing features and similarities.
    if (length(intersect(names(new.compound.features), rownames(existing.features))) == 0) {
        return(list(updated.features=existing.features, updated.similarities=existing.similarities))
    }
    
    # Add a new column to the feature matrix to hold the values of the features for 
    # the new compound.
    new.col <- rep(NA, nrow(existing.features))
    existing.features <- cbind(existing.features, new.col)
    colnames(existing.features)[ncol(existing.features)] <- new.compound.name
    
    # Now we can set the values of the new column according to new.compound.features.
    common.feature.names <- intersect(rownames(existing.features), names(new.compound.features))
    existing.features[common.feature.names, new.compound.name] <- new.compound.features[common.feature.names] 

    # Compute correlation of the new compound with existing compounds.
    new.similarities <- ComputePerturbationCorrelations(existing.features, new.compound.features,
                                                        common.feature.names)
    
    # Augment the similarity matrix to include the similarities with the new compound.
    existing.similarities <- rbind(existing.similarities, numeric(ncol(existing.similarities)))
    existing.similarities <- cbind(existing.similarities, numeric(nrow(existing.similarities)))
    
    colnames(existing.similarities)[ncol(existing.similarities)] <- new.compound.name
    rownames(existing.similarities)[nrow(existing.similarities)] <- new.compound.name
    
    existing.similarities[new.compound.name, ] <- new.similarities[colnames(existing.similarities)]
    existing.similarities[, new.compound.name] <- new.similarities[rownames(existing.similarities)]
    existing.similarities <- existing.similarities[order(rownames(existing.similarities)), 
                                                   order(colnames(existing.similarities))]
    
    list(updated.features=existing.features, updated.similarities=existing.similarities)
}

ComputePerturbationCorrelations <- function(existing.features, new.compound.features, common.feature.names) {
    # Compute the correlations between the features for the new compound and the existing features.
    #
    # Args:
    #   existing.features: A feature matrix where colnames correpond to drug names and rownames 
    #   correspond to gene names.
    #   new.compound.features: A named numeric vector containing features for the new compound.
    #   common.feature.names: A character vector containing features that are common between 
    #   the new compound and the existing features.
    #
    # Returns:
    #   A named numeric vector containing the similarities between the new compound and existing compounds.
    placeholder <- rep(NA, nrow(existing.features))
    names(placeholder) <- rownames(existing.features)
    
    placeholder[common.feature.names] <- new.compound.features[common.feature.names]
    new.similarities <- apply(existing.features, 2, cor, y = placeholder,
                              use="pairwise.complete.obs")
    
    new.similarities
}

LoadImagingFeatures <- function(file.path) {
    # Load imaging features from the specified file path.
    #
    # Args:
    #   file.path A character(1) containing the file path.
    #
    # Returns:
    #   Imaging features from the specified file path.
    imaging.features <- readRDS(file.path)
    
    imaging.features
}

LoadImagingSimilarities <- function(file.path) {
    # Load imaging similarities from the specified file path.
    #
    # Args:
    #   file.path A character(1) containing the file path.
    #
    # Returns:
    #   Imaging similarities from the specified file path.
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
    
    structure.features <- LoadStructureFeatures(feature.file.paths[['structure']])
    perturbation.features <- LoadPerturbationFeatures(feature.file.paths[['perturbation']])
    sensitivity.features <- LoadSensitivityFeatures(feature.file.paths[['sensitivity']])
    #imaging.features <- LoadImagingFeatures(feature.file.paths[['imaging']])

    existing.features <- list(sensitivity=sensitivity.features, perturbation=perturbation.features,
                              structure=structure.features, imaging=imaging.features)
    
    structure.similarities <- NULL
    perturbation.similarities <- NULL
    sensitivity.similarities <- NULL
    imaging.similarities <- NULL
        
    structure.similarities <- LoadStructureSimilarities(similarity.file.paths[['structure']])
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