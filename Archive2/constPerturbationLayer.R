###############################################################################################################
## Function reads in the perturbation data, and generates "affinity matrix" for the perturbation layer to be used either 
## solely or in combination of other layers (see function "integrateStrctSensPert") 
## input: 
##     sensDat: perturbation data processed in "perturbationData" function 
##     
## output: 
##     affinity matrix 
##
## 
###############################################################################################################


constPerturbationLayer <- function(pertDat) {

   # Correlation for Perturbation
   pertCor <- cor(pertDat, method = "pearson", use = "pairwise.complete.obs")
   ## Calculate affinity matrix (from generic distance) as described in SNFtool package with default values
   pertAff <- SNFtool::affinityMatrix(1-pertCor, 20, 0.5)

   return(pertAff)

}