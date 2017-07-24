###############################################################################################################
## Function reads in the perturbation data, processes/filters according to the intersection subset provided in the preprocess stage
##
## input: 
##     intersc: a vector of "character" (for "ctrpv2") or dataframe (for "nci60") containting intersection info.
## output: 
##     perturbation data ("dataframe") 	
##
## 
###############################################################################################################


PerturbationDataFlexibleCustomSig <- function(pert.file.name) {
    
    pert.data <- readRDS(pert.file.name)
    
    return(pert.data)
    
}