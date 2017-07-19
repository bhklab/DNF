CreateAugmentedMatrixSkeletons <- function(all.drugs) {
    augmented.matrices <- list(sens=NULL, pert=NULL, strc=NULL)
    
    for (i in 1:length(augmented.matrices)) {
        layer.name <- names(augmented.matrices)[i]
        
        augmented.matrices[[layer.name]] <- matrix(0, nrow=length(all.drugs), ncol=length(all.drugs))
        colnames(augmented.matrices[[layer.name]]) <- all.drugs
        rownames(augmented.matrices[[layer.name]]) <- all.drugs
    }
    
    augmented.matrices
}

ReplaceAugmentedExistingValues <- function(augmented.matrices, correlation.matrices) {
    for (k in 1:length(augmented.matrices)) {
        layer.name <- names(augmented.matrices)[k]
        
        drug.names <- rownames(correlation.matrices[[layer.name]])
        
        augmented.matrices[[layer.name]][drug.names, drug.names] <- correlation.matrices[[layer.name]]
    }

    augmented.matrices
}

CreateAffinityMatrices <- function(augmented.matrices) {
    affinity.matrices <- list(sens=NULL, pert=NULL, strc=NULL)
    
    for (i in 1:length(affinity.matrices)) {
        layer.name <- names(affinity.matrices)[i]
        
        affinity.matrices[[layer.name]] <- SNFtool::affinityMatrix(1 - augmented.matrices[[layer.name]], 20, 0.5)
    }
    
    affinity.matrices
}

ReplaceAffinityMatrixValues <- function(affinity.matrices, correlation.matrices, all.drugs) {
    for (i in 1:length(all.drugs)) {
        drug.name <- all.drugs[i]
        
        for (k in 1:length(affinity.matrices)) {
            layer.name <- names(affinity.matrices)[k]
            other.layers <- names(affinity.matrices)[-k]
            
            drugs.not.in.layer <- setdiff(all.drugs, colnames(correlation.matrices[[layer.name]]))
            
            if (drug.name %in% rownames(correlation.matrices[[layer.name]])) {
                for (d in drugs.not.in.layer) {
                    values.from.other.layers <- vapply(other.layers, function(x) {
                        if (d %in% rownames(correlation.matrices[[x]])) {
                            return(affinity.matrices[[x]][drug.name, d])
                        } else {
                            return(NA)
                        }
                        
                    }, numeric(1))                
                    
                    
                    values.from.other.layers <- na.omit(values.from.other.layers)
                    
                    imputed.value <- mean(values.from.other.layers)
                    
                    affinity.matrices[[layer.name]][drug.name, d] <- imputed.value
                    affinity.matrices[[layer.name]][d, drug.name] <- imputed.value
                }            
            } else {
                for (d in all.drugs) {
                    values.from.other.layers <- vapply(other.layers, function(x) {
                        if (d %in% rownames(correlation.matrices[[x]])) {
                            return(affinity.matrices[[x]][drug.name, d])
                        } else {
                            return(NA)
                        }
                        
                    }, numeric(1))                
                    
                    
                    values.from.other.layers <- na.omit(values.from.other.layers)
                    
                    imputed.value <- mean(values.from.other.layers)
                    
                    affinity.matrices[[layer.name]][drug.name, d] <- imputed.value
                    affinity.matrices[[layer.name]][d, drug.name] <- imputed.value
                    
                }
            }
        }
    }
    
    affinity.matrices
}