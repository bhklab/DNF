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


perturbationData <- function(dname="lincs", intersc) {

    ## Keep drugs from LINCS intersecting with ctrpv2 
    if (class(intersc)=="character") {
        lincs <- read.csv("Data/LINCS.csv",stringsAsFactors = FALSE) # 20326 drugs x 28 descriptions
        ## capitalize + remobe badchars and intersect LINCS chemicals from LINCS with ctrpv2 compounds
        lincs$pert_iname <- toupper(lincs$pert_iname)
        lincs$pert_iname <- gsub(badchars,"",lincs$pert_iname)
        lincsInters <- lincs[lincs$pert_iname %in% intersc,,drop=F]
        lincsInters <- lincsInters[!duplicated(lincsInters$pert_iname),,drop=F] ## 239 drugs 
     } else { ## class is dataframe ... for nci60
        lincsInters <- intersc
     }
  
     ## read the perturbation file
     load("Data/l1000_drug_signatures.RData")
     pertLincs <- l1000.drug.signatures
     ## extract estimates of the drug pert signatures
     dfLincs <- pertLincs[,,"estimate"] # 978 genes x 20364 drugs
     ## load pheno data and match with gene profiles
     dfLincs <- dfLincs[, match(lincsInters$pert_id,colnames(dfLincs))] #239 drugs
     ## genes symbol replacement from genids
     rownames(dfLincs) <- gsub(x = rownames(dfLincs),pattern = "geneid.",replacement = "")
     symb <- annotate::lookUp(rownames(dfLincs), 'org.Hs.eg', 'SYMBOL')
     rownames(dfLincs) <- unlist(unname(symb)) 
     dfLincs <- dfLincs[!is.na(rownames(dfLincs)),]
     if (!(all(colnames(dfLincs) == lincsInters$pert_id))) { #Sanity check
        stop(paste("error!")) }
     colnames(dfLincs) <- lincsInters$pert_iname
     dfLincs <- dfLincs[, !apply(is.na(dfLincs), 2, all)] #Remove drugs which has all NA in columns
     dim(dfLincs) #237 drugs after all filtering
     dfLincs <- dfLincs[, order(colnames(dfLincs)),drop=F]

     return(dfLincs)

}