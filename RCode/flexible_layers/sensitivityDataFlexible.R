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



SensitivityDataFlexible <- function(combined_file_path="") {
    combinedSens <- readRDS(combined_file_path)
    combinedSens <- combinedSens[!duplicated(rownames(combinedSens)), , drop=FALSE]
    sens <- combinedSens

    return(sens)
    
}