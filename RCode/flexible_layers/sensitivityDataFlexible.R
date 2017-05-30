###############################################################################################################
## Function reads in the sensitivity data, processes/filters according to the intersection subset provided in the preprocess stage
##
## input: 
##     dname: name ("character") of the dataset, e.g., "nci60" , "ctrpv2"   
##     intersc: a vector of "character" (for "ctrpv2") or dataframe (for "nci60") containting intersection info.
## output: 
##     sensitivity data ("dataframe") 	
##
## 
###############################################################################################################



SensitivityDataFlexible <- function(combined.file.path="") {
    combined.sens <- readRDS(combined.file.path)
    combined.sens <- combined.sens[!duplicated(rownames(combined.sens)), , drop=FALSE]
    sens <- combined.sens

    return(sens)
    
}