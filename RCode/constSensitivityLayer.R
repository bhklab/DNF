###############################################################################################################
## Function reads in the sensitivity data, and generates "affinity matrix" for the sensitivity layer to be used either 
## solely or in combination of other layers (see function "integrateStrctSensPert") 
## input: 
##     sensDat: sesitivity data processed in "sensitivityData" function 
##     
## output: 
##     affinity matrix 
##
## 
###############################################################################################################


constSensitivityLayer <- function(sensDat) {
  
  # Correlation for Sensivity
  sensCor <- cor(sensDat, method = "pearson", use = "pairwise.complete.obs")
  ## if NA remaining in cor matrix, replace with 0s, not very clean but no other choices for now
  sensCor <- apply(sensCor, 1, function(x) ifelse(is.na(x),0,x))
  ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
  sensAff <- SNFtool::affinityMatrix(1-sensCor, 20, 0.5)
  
  return(sensAff)
  
}