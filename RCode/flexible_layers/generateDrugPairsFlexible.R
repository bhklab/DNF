###############################################################################################################
## Function intersects the benchmark data and the layers of the network and then generates drug pairs for the 
## final (intersection) subset of drugs
##
## input: 
##     benchDat
##     benchDat
##     strcAff
##     sensAff
##     pertAff
##     integration  
##
## output: 
##     list of drug pairs from benchmark, structure, perturbation, sensitivtiy and the integration of all
##
## 
###############################################################################################################

GenerateDrugPairsFlexible <- function(benchDat, strcAff=NULL, sensAff=NULL, 
                                    pertAff=NULL, integration=NULL, 
                                    luminexAff=NULL, imagingAff=NULL) {
    
    layers <- list(strcAff=strcAff, sensAff=sensAff, pertAff=pertAff, integr=integration,
                   luminexAff=luminexAff, imagingAff=imagingAff)
    
    layers <- layers[!sapply(layers, is.null)]
    intx <- intersect(colnames(benchDat), colnames(layers[[1]]))
    
    data.list <- list()
    
    for (i in 1:length(layers)) {
        prefix <- strsplit(names(layers)[[i]], split="[[:upper:]]")[[1]][1]
        prefix <- paste0(toupper(substr(prefix, 1, 1)), substr(prefix, 2, nchar(prefix)))
        
        dataset.name <- paste("data", prefix, sep="")
        aff.mat <- layers[[i]]
        
        temp <- aff.mat[intx, intx] 
        temp[upper.tri(temp, diag=TRUE)] <- NA
        data.list[[dataset.name]] <- temp 
    }
    
    pairs.list <- list()
    
    benchDat[upper.tri(benchDat, diag=TRUE)] <- NA
    benchPairs <- melt(benchDat)
    benchPairs <- na.omit(benchPairs)
    colnames(benchPairs)[3] <- "bench"
    
    for (i in 1:length(data.list)) {
        prefix <- strsplit(names(data.list)[i], split="data")[[1]][2]
        prefix <- tolower(prefix)
        
        pairs.name <- paste(prefix, "Pairs", sep="")
        temp <- melt(data.list[[i]])
        temp <- na.omit(temp)
        colnames(temp)[3] <- paste("obs", prefix, sep=".")
        
        pairs.list[[pairs.name]] <-  temp
    }
    
    pairs.list[["benchPairs"]] <- benchPairs
    
    return(pairs.list)
}