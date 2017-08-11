IntegrateImputedLayers <- function(sens.cor, pert.cor, strc.cor, all.drugs) {
    correlation.matrices <- list(sens=sens.cor, pert=pert.cor, strc=strc.cor)
    
    affinity.matrices <- CreateAffinityMatrices(correlation.matrices)
    augmented.matrices <- CreateAugmentedMatrixSkeletons(names(correlation.matrices), all.drugs)
    augmented.matrices <- ReplaceAugmentedExistingValues(augmented.matrices, affinity.matrices)
    affinity.matrices <- ReplaceAffinityMatrixValuesFast(augmented.matrices, correlation.matrices, 
                                                         all.drugs)
    affinity.matrices <- medianSimilarity(affinity.matrices)
    
    integrated <- SNFtool::SNF(affinity.matrices)
    rownames(integrated) <- all.drugs
    colnames(integrated) <- all.drugs
    
    return(list(integrated=integrated, affinity.matrices=affinity.matrices))
}