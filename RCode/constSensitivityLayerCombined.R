###############################################################################################################
## Function takes sensitivity data, and generates "affinity matrix" for the sensitivity layer to be used either 
## solely or in combination of other layers (see function "integrateStrctSensPert") 
## input: 
##     sensMat: sesitivity matrix obtained in "sensitivityData" function with "combined" argument
##     
## output: 
##     affinity matrix 
##
## 
###############################################################################################################


constSensitivityLayerCombined <- function(sensMat) {
    ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
    sensAff <- SNFtool::affinityMatrix(1-sensMat, 20, 0.5)
    
    return(sensAff)
    
}