###############################################################################################################
## Function takes sensitivity data, and generates "affinity matrix" for the sensitivity layer to be used either 
## solely or in combination of other layers (see function "integrateStrctSensPert") 
## input: 
##     sens.mat: sesitivity matrix obtained in "sensitivityData" function with "combined" argument
##     
## output: 
##     affinity matrix 
##
## 
###############################################################################################################


ConstSensitivityLayerFlexible <- function(sens.mat) {
    ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
    sens.aff <- SNFtool::affinityMatrix(1-sens.mat, 20, 0.5)
    
    return(sens.aff)
    
}