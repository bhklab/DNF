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


integrateStrctSensPert <- function(sensAff, strcAff, pertAff) {
  
  integration <- SNFtool::SNF(list(sensAff, strcAff, pertAff))
  colnames(integration) <- rownames(integration) <- colnames(strcAff)
  
  return(integration)
  
}

