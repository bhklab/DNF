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


PerturbationDataFlexible <- function(pert.file.name, lincs.meta, use.subsetted=FALSE) {
    if (use.subsetted) {
        df.lincs <- readRDS(pert.file.name)
        return(df.lincs)
    }
    
    
    ## read the perturbation file
    load(pert.file.name)
    pert.lincs <- L1000_compounds.perturbation        

    
    ## extract estimates of the drug pert signatures
    df.lincs <- pert.lincs[,,"estimate"] # 978 genes x 20364 drugs. Actually looks like it's 978 genes x 414 drugs
    ## load pheno data and match with gene profiles
    df.lincs <- df.lincs[, match(lincs.meta$pert_id,colnames(df.lincs))] #239 drugs
    ## genes symbol replacement from genids
    rownames(df.lincs) <- gsub(x = rownames(df.lincs),pattern = "geneid.",replacement = "")
    symb <- annotate::lookUp(rownames(df.lincs), 'org.Hs.eg', 'SYMBOL')
    rownames(df.lincs) <- unlist(unname(symb)) 
    df.lincs <- df.lincs[!is.na(rownames(df.lincs)),]
    
    colnames(df.lincs) <- lincs.meta$pert_iname
    df.lincs <- df.lincs[, !apply(is.na(df.lincs), 2, all)] #Remove drugs which has all NA in columns
    dim(df.lincs) #237 drugs after all filtering. ACtually looks like it's 239
    df.lincs <- df.lincs[, order(colnames(df.lincs)),drop=F]
    
    return(df.lincs)
    
}