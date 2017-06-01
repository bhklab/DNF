AverageAUCS <- function(overlap, aucs.all, aucs.all.final) {
    for (i in 1:length(overlap)) {
        drug <- overlap[i]
        
        drug.aucs <- aucs.all[endsWith(rownames(aucs.all), drug), ]
        if (class(drug.aucs) == "data.frame" | class(drug.aucs) == "matrix") {
            combined.aucs <- colMeans(drug.aucs, na.rm = TRUE)    
        } else {
            combined.aucs <- drug.aucs
        }
        
        #temp <- matrix(combined.aucs, nrow=1, ncol=length(drug.aucs))
        #temp <- as.data.frame(temp)
        
        aucs.all.final[i, ] <- combined.aucs
    } 
    
    aucs.all.final
}

