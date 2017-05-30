###############################################################################################################
## Function reads in the affinity matrices or network layers of sensitivity, structure and perturbation data, 
## and integrates them via SNF (similarity network fusion) method 
## input: 
##     sens.aff: affinity matrix generated in the "constSensitivityLayer"
##     strc.aff: affinity matrix generated in the "constStructureLayer"
##     pert.aff: affinity matrix generated in the "constPerturbationLayer"
##     
## output: 
##     integartion result (matrix)  
##
## 
###############################################################################################################


IntegrateLayersFlexible <- function(sens.aff=NULL, strc.aff=NULL, pert.aff=NULL,
                                    luminex.aff=NULL, imaging.aff=NULL) {
    layers <- list(sens.aff=sens.aff, strc.aff=strc.aff, pert.aff=pert.aff,
                   luminex.aff=luminex.aff, imaging.aff=imaging.aff)
    layers <- layers[!sapply(layers, is.null)]
    
    if (length(layers) > 1) {
        integration <- SNFtool::SNF(layers)   
        colnames(integration) <- rownames(integration) <- colnames(layers[[1]])
    } else {
        integration <- layers[[1]]
    }
    

    return(integration)
}

