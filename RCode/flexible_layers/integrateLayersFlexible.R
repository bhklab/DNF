###############################################################################################################
## Function reads in the affinity matrices or network layers of sensitivity, structure and perturbation data, 
## and integrates them via SNF (similarity network fusion) method 
## input: 
##     sensAff: affinity matrix generated in the "constSensitivityLayer"
##     strcAff: affinity matrix generated in the "constStructureLayer"
##     pertAff: affinity matrix generated in the "constPerturbationLayer"
##     
## output: 
##     integartion result (matrix)  
##
## 
###############################################################################################################


IntegrateLayersFlexible <- function(sensAff=NULL, strcAff=NULL, pertAff=NULL, luminexAff=NULL, imagingAff=NULL) {
    layers <- list(sensAff=sensAff, strcAff=strcAff, pertAff=pertAff,
                   luminexAff=luminexAff, imagingAff=imagingAff)
    layers <- layers[!sapply(layers, is.null)]
    
    if (length(layers) > 1) {
        integration <- SNFtool::SNF(layers)   
        colnames(integration) <- rownames(integration) <- colnames(layers[[1]])
    } else {
        integration <- layers[[1]]
    }
    

    return(integration)
    
}

