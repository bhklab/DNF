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


integrateStrctSensPertModded <- function(sensAff, strcAff, pertAff, luminexAff=NULL, imagingAff=NULL) {
    if (is.null(luminexAff) && is.null(imagingAff)) {
        integration <- SNFtool::SNF(list(sensAff, strcAff, pertAff))    
    } else if (is.null(imagingAff)) {
        integration <- SNFtool::SNF(list(sensAff, strcAff, pertAff, luminexAff))    
    } else if (is.null(luminexAff)) {
        integration <- SNFtool::SNF(list(sensAff, strcAff, pertAff, imagingAff))    
    } else {
        integration <- SNFtool::SNF(list(sensAff, strcAff, pertAff, imagingAff, luminexAff))    
    }
    
    colnames(integration) <- rownames(integration) <- colnames(strcAff)
    
    return(integration)
    
}

