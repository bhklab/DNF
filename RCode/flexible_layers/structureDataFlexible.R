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


StructureDataFlexible <- function(lincs.meta.subset) 
{
    ## Matrix of SMILES (tanimoto distance) 
    ## parse smiles 
    targetMolcs <- rcdk::parse.smiles(as.character(lincs.meta.subset$canonical_smiles))  
    ## assign fingerprints (extended connectivity fingerprint)
    targetFps <- lapply(targetMolcs, rcdk::get.fingerprint, type='extended') 
    names(targetFps) <- lincs.meta.subset$pert_iname
    length(targetFps) 
    res <- targetFps
    
    return(res)
}