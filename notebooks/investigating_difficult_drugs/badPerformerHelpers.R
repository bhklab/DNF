IntegrateCorrelationMatrices <- function(correlation.matrices, all.drugs) {
    affinity.matrices <- CreateAffinityMatrices(correlation.matrices)
    augmented.matrices <- CreateAugmentedMatrixSkeletons(names(correlation.matrices), all.drugs)
    augmented.matrices <- ReplaceAugmentedExistingValues(augmented.matrices, affinity.matrices)
    affinity.matrices <- ReplaceAffinityMatrixValuesFast(augmented.matrices, correlation.matrices,
                                                         all.drugs)
    affinity.matrices <- medianSimilarity(affinity.matrices)
    
    integrated <- SNFtool::SNF(affinity.matrices)
    rownames(integrated) <- all.drugs
    colnames(integrated) <- all.drugs
    
    return(integrated)
}