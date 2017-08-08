GetImputedAccuracies <- function(to.impute, drugs.to.remove, pert.cor, sens.cor, strc.cor) {
    c1 <- makeCluster(8)
    registerDoParallel(c1)
    
    k <- 7
    accuracies <- list()
    accuracies <- foreach(i = names(drugs.to.remove), .final = function(x) setNames(x, names(drugs.to.remove))) %do%
    {
        require(foreach)
        accuracies[[as.character(i)]] <- c()
        
        foreach(j = 1:length(drugs.to.remove[[i]])) %dopar% {
            to.remove <- drugs.to.remove[[as.character(i)]][[j]]
            
            sub.sens.cor <- sens.cor
            sub.pert.cor <- pert.cor
            sub.strc.cor <- strc.cor
            
            if (to.impute == 'sens') {
                sub.sens.cor <- sens.cor[!(rownames(sens.cor) %in% to.remove), !(colnames(sens.cor) %in% to.remove)]    
            } else if (to.impute == 'pert') {
                sub.pert.cor <- pert.cor[!(rownames(pert.cor) %in% to.remove), !(colnames(pert.cor) %in% to.remove)]
            } else if (to.impute == 'strc') {
                sub.strc.cor <- strc.cor[!(rownames(strc.cor) %in% to.remove), !(colnames(strc.cor) %in% to.remove)]
            }
            
            correlation.matrices <- list(sens=sub.sens.cor, pert=sub.pert.cor, strc=sub.strc.cor)
            
            affinity.matrices <- CreateAffinityMatrices(correlation.matrices)
            augmented.matrices <- CreateAugmentedMatrixSkeletons(names(correlation.matrices), all.drugs)
            augmented.matrices <- ReplaceAugmentedExistingValues(augmented.matrices, affinity.matrices)
            affinity.matrices <- ReplaceAffinityMatrixValuesFast(augmented.matrices, correlation.matrices,
                                                                 all.drugs)
            affinity.matrices <- medianSimilarity(affinity.matrices)
            
            integrated <- SNFtool::SNF(affinity.matrices)
            rownames(integrated) <- all.drugs
            colnames(integrated) <- all.drugs
            
            accuracy <- GetKNNAccuracy(k, integrated, data.bench)
            
            accuracy <- sum(accuracy) / length(accuracy)
            
            accuracy
        }
    }
    
    stopCluster(c1)
    return(accuracies)
}

CreateDataForGGPlot <- function(acccuracies.pert, accuracies.sens, accuracies.strc) {
    accuracies.unlisted.sens <- lapply(accuracies.sens, unlist)
    accuracies.unlisted.pert <- lapply(accuracies.pert, unlist)
    accuracies.unlisted.strc <- lapply(accuracies.strc, unlist)
    
    accuracies.means.sens <- unlist(lapply(accuracies.unlisted.sens, mean))
    accuracies.means.pert <- unlist(lapply(accuracies.unlisted.pert, mean))
    accuracies.means.strc <- unlist(lapply(accuracies.unlisted.strc, mean))
    
    accuracies.std.sens <- unlist(lapply(accuracies.unlisted.sens, sd))
    accuracies.std.strc <- unlist(lapply(accuracies.unlisted.strc, sd))
    accuracies.std.pert <- unlist(lapply(accuracies.unlisted.pert, sd))
    accuracies.std.sens[is.na(accuracies.std.sens)] <- 0
    accuracies.std.strc[is.na(accuracies.std.strc)] <- 0
    accuracies.std.pert[is.na(accuracies.std.pert)] <- 0
    
    
    accuracy.stats <- data.frame(sens.mean=accuracies.means.sens, sens.std=accuracies.std.sens,
                                 strc.mean=accuracies.means.strc, strc.std=accuracies.std.strc,
                                 pert.mean=accuracies.means.pert, pert.std=accuracies.std.pert, n=as.numeric(names(accuracies.means)))   
    
    return(accuracy.stats)
}
