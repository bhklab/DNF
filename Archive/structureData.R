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
    if (class(intersc) == "character") { 
       lincs <- read.csv("Data/LINCS.csv",stringsAsFactors = FALSE) # 20326 drugs x 28 descriptions
       ## capitalize + remobe badchars and intersect LINCS chemicals from LINCS with ctrpv2 compounds
       lincs$pert_iname <- toupper(lincs$pert_iname)
       lincs$pert_iname <- gsub(badchars,"",lincs$pert_iname)
       lincsInters <- lincs[lincs$pert_iname %in% intersc,,drop=F]
       lincsInters <- lincsInters[!duplicated(lincsInters$pert_iname),,drop=F] ## 239 drugs 
     } else { 
       lincsInters <- intersc
     }
  
     ## Matrix of SMILES (tanimoto distance) 
     ## parse smiles 
     targetMolcs <- rcdk::parse.smiles(as.character(lincsInters$canonical_smiles))  
     ## assign fingerprints (extended connectivity fingerprint)
     targetFps <- lapply(targetMolcs, rcdk::get.fingerprint, type='extended') 
     names(targetFps) <- lincsInters$pert_iname
     length(targetFps) 
     res <- targetFps
     return(res)
     
  } else {
    stop(paste("please check the first argument."))
  }
  
  
}