###############################################################################################################
## Function reads in the structure data, processes/filters according to the intersection subset provided in the preprocess stage
##
## input: 
##     dname: currently must be set to "lincs"
##     intersc: a vector of "character" or dataframe containting intersection info.
## output: 
##     structure data ("dataframe") 	
##
## 
###############################################################################################################


structureData <- function(dname="lincs", intersc) 
{
  
  if (dname == "lincs") {
     ## Matrix of SMILES (tanimoto distance) 
     ## parse smiles 
     targetMolcs <- rcdk::parse.smiles(as.character(intersc$canonical_smiles))  
     ## assign fingerprints (extended connectivity fingerprint)
     targetFps <- lapply(targetMolcs, rcdk::get.fingerprint, type='extended') 
     names(targetFps) <- intersc$pert_iname
     length(targetFps) 
     res <- targetFps
     return(res)
     
  } else {
    stop(paste("please check the first argument."))
  }
  
  
}