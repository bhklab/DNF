###############################################################################################################
## Function reads in the perturbation data, processes/filters according to the intersection subset provided in the preprocess stage
##
## input: 
##     dname: currently must be set to "lincs"  
##     intersc: a vector of "character" (for "ctrpv2") or dataframe (for "nci60") containting intersection info.
## output: 
##     perturbation data ("dataframe") 	
##
## 
###############################################################################################################


perturbationData <- function(dname="lincs", intersc, name) {
  
     ## read the perturbation file
     load("Data/DNF_pert_all.RData")
     pertLincs <- l1000.drug.signatures
  
      ## extract estimates of the drug pert signatures
     dfLincs <- pertLincs[,,"estimate"] # 978 genes x 20364 drugs
     ## load pheno data and match with gene profiles
     dfLincs <- dfLincs[, match(intersc$pert_id,colnames(dfLincs))] #239 drugs
     ## genes symbol replacement from genids
     rownames(dfLincs) <- gsub(x = rownames(dfLincs),pattern = "geneid.",replacement = "")
     symb <- annotate::lookUp(rownames(dfLincs), 'org.Hs.eg', 'SYMBOL')
     rownames(dfLincs) <- unlist(unname(symb)) 
     dfLincs <- dfLincs[!is.na(rownames(dfLincs)),]
     #if (!(all(colnames(dfLincs) == intersc$pert_id))) { #Sanity check
      #  stop(paste("error!")) }
     colnames(dfLincs) <- intersc$pert_iname
     dfLincs <- dfLincs[, !apply(is.na(dfLincs), 2, all)] #Remove drugs which has all NA in columns
     dim(dfLincs) #237 drugs after all filtering
     dfLincs <- dfLincs[, order(colnames(dfLincs)),drop=F]

     return(dfLincs)

}