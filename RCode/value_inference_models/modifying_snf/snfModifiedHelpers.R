CreateAugmentedMatrixSkeletons <- function(matrix.names, all.drugs) {
    augmented.matrices <- list()
    
    for (i in 1:length(matrix.names)) {
        layer.name <- matrix.names[i]
        
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
    affinity.matrices <- list()
    
    for (i in 1:length(augmented.matrices)) {
        layer.name <- names(augmented.matrices)[i]
        
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
                    
                    if (!is.na(imputed.value)) {
                        affinity.matrices[[layer.name]][drug.name, d] <- imputed.value
                        affinity.matrices[[layer.name]][d, drug.name] <- imputed.value   
                    }
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
                    
                    if (!is.na(imputed.value)) {
                        affinity.matrices[[layer.name]][drug.name, d] <- imputed.value
                        affinity.matrices[[layer.name]][d, drug.name] <- imputed.value   
                    }                    
                }
            }
        }
    }
    
    affinity.matrices
}

ReplaceAffinityMatrixValuesFast <- function(affinity.matrices, correlation.matrices, all.drugs) {
    for (i in 1:length(all.drugs)) {
        drug.name <- all.drugs[i]
        
        for (k in 1:length(affinity.matrices)) {
            layer.name <- names(affinity.matrices)[k]
            other.layers <- names(affinity.matrices)[-k]
            
            drugs.not.in.layer <- setdiff(all.drugs, colnames(correlation.matrices[[layer.name]]))
 
            if (drug.name %in% rownames(correlation.matrices[[layer.name]])) {
                drugs.to.replace <- drugs.not.in.layer
            } else {
                drugs.to.replace <- all.drugs
            }
            
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