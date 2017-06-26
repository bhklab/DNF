CompConcordIndxFlexible <- function(all.pairs, compute.p.val=TRUE) {
    # Computes concordance indices and p-values for the integrative method vs. single layers.
    #
    # Args:
    #   all.pairs: A list containing dataframes where each dataframe has 3 columns: Var1,
    #              Var2, and obs.integr Var1 and Var2 are drug names, and obs.integr which is the 
    #              integrated similarity measure for the two drugs. One of the dataframes is 
    #              based on the drug target data, and instead of a obs.integr column, it has 
    #              a bench column which takes values in {0,1} indicating whether or not drugs
    #              share the same target.
    #   comp.p.val: A boolean indicating whether or not to compute p-values.
    #
    # Returns:
    #   A list containing 2 elements. The first element is itself a list of concordance indices
    #   (one for every layer on the model). The second element is list of p-values indicating
    #   whether the significance in improvement of the integration of layers vs. single layers.
    # 
    c.index.obj.list <- list()
    
    # Create a list of c index objects
    for (i in 1:length(all.pairs)) {
        pair.name <- strsplit(names(all.pairs)[i], "Pairs")[[1]][1]
        
        if (pair.name != "bench") {
            c.index.name <- paste(pair.name, "Cindex", sep="")
            
            temp <- Hmisc::rcorr.cens(x=all.pairs[[i]][, paste("obs", pair.name, sep=".")], S=all.pairs$benchPairs$bench)
            c.index.obj.list[[c.index.name]] <- temp
        }
    
    }
    c.index.list <- list()
    
    # Use the c.index.obj.list to obtain the actual C index values from the objects
    for (i in 1:length(c.index.obj.list)) {
        c.index.list[[names(c.index.obj.list)[i]]] <- c.index.obj.list[[i]]["C Index"]
    }
    
    p.vals.list <- list()
    
    # Compute the p-values to see if the improvement of the integration over the single layers
    # is significant.
    if (compute.p.val) {
        for (i in 1:length(c.index.obj.list)) {
            c.index.name <- names(c.index.obj.list)[i]
            
            if (c.index.name != "integrCindex") {
                prefix <- strsplit(c.index.name, "Cindex")[[1]][1]
                
                pair.name <- paste(prefix, "Pairs", sep="")
                col.name <- paste("obs", prefix, sep=".")
                
                temp <- cindexComp2(c.index.obj.list$integrCindex, c.index.obj.list[[c.index.name]],
                                    all.pairs$integrPairs$obs.integr, all.pairs[[pair.name]][col.name])
                temp <- temp$p.value
                
                p.val.name <- paste(prefix, "PVal", sep="")
                p.vals.list[[p.val.name]] <- temp
            }
        }        
    }
    
    # Return both lists of results
    r <- list(c.index.list=c.index.list, p.vals.list=p.vals.list) 
    
    return(r)
}


