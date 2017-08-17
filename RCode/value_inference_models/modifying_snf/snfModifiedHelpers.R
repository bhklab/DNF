CreateAugmentedMatrixSkeletons <- function(matrix.names, all.drugs) {
    augmented.matrices <- list()
    
    for (i in 1:length(matrix.names)) {
        layer.name <- matrix.names[i]
        
        augmented.matrices[[layer.name]] <- matrix(NA, nrow=length(all.drugs), ncol=length(all.drugs))
        colnames(augmented.matrices[[layer.name]]) <- all.drugs
        rownames(augmented.matrices[[layer.name]]) <- all.drugs
    }
    
    augmented.matrices
}

ReplaceAugmentedExistingValues <- function(augmented.matrices, correlation.matrices) {
    for (k in 1:length(augmented.matrices)) {
        layer.name <- names(augmented.matrices)[k]
        
        drug.names <- rownames(correlation.matrices[[layer.name]])
        drug.names <- intersect(rownames(augmented.matrices[[layer.name]]), drug.names)
        
        # correlation.matrices[[layer.name]][is.na(correlation.matrices[[layer.name]])] <- 0
        
        augmented.matrices[[layer.name]][drug.names, drug.names] <- correlation.matrices[[layer.name]][drug.names, drug.names]
    }

    augmented.matrices
}

CreateAffinityMatrices <- function(augmented.matrices) {
    affinity.matrices <- list()
    
    for (i in 1:length(augmented.matrices)) {
        layer.name <- names(augmented.matrices)[i]
        
        affinity.matrices[[layer.name]] <- SNFtool::affinityMatrix(1 - augmented.matrices[[layer.name]], 20, 0.5)
    }
    
    affinity.matrices
}

ReplaceAffinityMatrixValuesFast <- function(affinity.matrices, correlation.matrices, all.drugs) {
    if (length(affinity.matrices) == 1) {
        return(affinity.matrices)
    }
    
    for (i in 1:length(all.drugs)) {
        drug.name <- all.drugs[i]
        
        for (k in 1:length(affinity.matrices)) {
            layer.name <- names(affinity.matrices)[k]
            other.layers <- names(affinity.matrices)[-k]
            
            drugs.not.in.layer <- setdiff(all.drugs, colnames(correlation.matrices[[layer.name]]))
 
                
            if (drug.name %in% rownames(correlation.matrices[[layer.name]])) {
                # When this drug is present in the layer, all we have to do is overwrite
                # the affinity with this drug and the drugs that are not in the layer
                drugs.to.replace <- drugs.not.in.layer
            } else {
                # When this drug is missing in the present layer, we have to overwrite the affinity 
                # with this drug and every other possible drug
                drugs.to.replace <- all.drugs
            }
            
            # For the layers that the current drug can be found in, determine which drugs
            # from each layer can be used in overwriting the
            in.layers <- lapply(correlation.matrices[-k], function(x) {
                if (drug.name %in% rownames(x)) {
                    intersect(drugs.to.replace, rownames(x))    
                } else {
                    character(0)
                }
            })
            
            values.to.replace <- lapply(seq_along(in.layers), function(x, y, layer.names) {
                temp <- rep(NA, length(drugs.to.replace))
                names(temp) <- drugs.to.replace
                
                temp.layer.name <- layer.names[x]
                names.of.interest <- y[[x]]
                
                temp[names.of.interest] <- affinity.matrices[[temp.layer.name]][drug.name, names.of.interest]
                
                temp
                
            }, y=in.layers, layer.names = names(in.layers))
            
            values.to.replace <- do.call(rbind, values.to.replace)
            values.to.replace <- colMeans(values.to.replace, na.rm=TRUE)
            names(values.to.replace) <- drugs.to.replace
            
            nan.names <- names(values.to.replace[is.nan(values.to.replace)])
            values.to.replace[nan.names] <- affinity.matrices[[layer.name]][drug.name, nan.names]
            
            affinity.matrices[[layer.name]][drug.name, drugs.to.replace] <- values.to.replace[drugs.to.replace]
            affinity.matrices[[layer.name]][drugs.to.replace, drug.name] <- values.to.replace[drugs.to.replace]
                        
        }
    }
    
    affinity.matrices
}

ReplaceCorrelationValues <- function(augmented.matrices, correlation.matrices, all.drugs) {
    for (i in 1:length(all.drugs)) {
        drug.name <- all.drugs[i]
        
        for (k in 1:length(augmented.matrices)) {
            layer.name <- names(augmented.matrices)[k]
            other.layers <- names(augmented.matrices)[-k]
            
            drugs.not.in.layer <- setdiff(all.drugs, colnames(correlation.matrices[[layer.name]]))
            
            if (drug.name %in% rownames(correlation.matrices[[layer.name]])) {
                # When this drug is present in the layer, all we have to do is overwrite
                # the affinity with this drug and the drugs that are not in the layer
                drugs.to.replace <- drugs.not.in.layer
            } else {
                # When this drug is missing in the present layer, we have to overwrite the affinity 
                # with this drug and every other possible drug
                drugs.to.replace <- all.drugs
            }
            
            # For the layers that the current drug can be found in, determine which drugs
            # from each layer can be used in overwriting the
            in.layers <- lapply(correlation.matrices[-k], function(x) {
                if (drug.name %in% rownames(x)) {
                    intersect(drugs.to.replace, rownames(x))    
                } else {
                    character(0)
                }
            })
            
            values.to.replace <- lapply(seq_along(in.layers), function(x, y, layer.names) {
                temp <- rep(NA, length(drugs.to.replace))
                names(temp) <- drugs.to.replace
                
                temp.layer.name <- layer.names[x]
                names.of.interest <- y[[x]]
                
                temp[names.of.interest] <- augmented.matrices[[temp.layer.name]][drug.name, names.of.interest]
                
                temp
                
            }, y=in.layers, layer.names = names(in.layers))
            
            values.to.replace <- do.call(rbind, values.to.replace)
            values.to.replace <- colMeans(values.to.replace, na.rm=TRUE)
            names(values.to.replace) <- drugs.to.replace
            
            # This part is where we end up saying that for the drugs that we couldn't get data
            # from other layers with which to overwrite the current layer, simply keep the values of
            # the current layer. In this case, that means 0. This might not be the best thing to do since
            # we're saying for sure that the current drug and the drugs in nan.names have a similarity of 0.
            nan.names <- names(values.to.replace[is.nan(values.to.replace)])
            values.to.replace[nan.names] <- augmented.matrices[[layer.name]][drug.name, nan.names]
            
            augmented.matrices[[layer.name]][drug.name, drugs.to.replace] <- values.to.replace[drugs.to.replace]
            augmented.matrices[[layer.name]][drugs.to.replace, drug.name] <- values.to.replace[drugs.to.replace]
        }
    }
    
    augmented.matrices
}

IntegrateCorrelationMatrices <- function(correlation.matrices, all.drugs) {
    affinity.matrices <- CreateAffinityMatrices(correlation.matrices)
    augmented.matrices <- CreateAugmentedMatrixSkeletons(names(correlation.matrices), all.drugs)
    augmented.matrices <- ReplaceAugmentedExistingValues(augmented.matrices, affinity.matrices)
    affinity.matrices <- ReplaceAffinityMatrixValuesFast(augmented.matrices, correlation.matrices,
                                                         all.drugs)
    
    # Any remaining NAs simply get replaced via median value imputation
    affinity.matrices <- medianSimilarity(affinity.matrices)
    
    if (length(affinity.matrices) > 1) {
        integrated <- SNFtool::SNF(affinity.matrices)
        rownames(integrated) <- all.drugs
        colnames(integrated) <- all.drugs
    } else {
        integrated <- affinity.matrices[[1]]
    }

    return(integrated)
}