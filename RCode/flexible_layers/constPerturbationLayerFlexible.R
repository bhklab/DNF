###############################################################################################################
## Function reads in the perturbation data, and generates "affinity matrix" for the perturbation layer to be used either 
## solely or in combination of other layers (see function "integrateStrctSensPert") 
## input: 
##     pert.data: perturbation data processed in "perturbationData" function 
##     
## output: 
##     affinity matrix 
##
## 
###############################################################################################################


ConstPerturbationLayerFlexible <- function(pert.data) {
    # Correlation for Perturbation
    pert.cor <- cor(pert.data, method = "pearson", use = "pairwise.complete.obs")
    pert.cor[is.na(pert.cor)] <- 0
    
    saveRDS(pert.cor, "Data/upload/similarities/perturbation/perturbation_similarities.RData")
    ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
    pert.aff <- SNFtool::affinityMatrix(1-pert.cor, 20, 0.5)
    
    return(pert.aff)
    
}