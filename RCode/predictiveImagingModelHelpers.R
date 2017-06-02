CreateAUCFinalMatrix <- function(aucs.all, overlap) {
    aucs.all.final <- matrix(0, ncol=ncol(aucs.all), nrow=length(overlap))
    rownames(aucs.all.final) <- overlap
    colnames(aucs.all.final) <- colnames(aucs.all)
    
    aucs.all.final <- AverageAUCS(overlap, aucs.all, aucs.all.final)
    columns.to.keep <- colSums(is.na(aucs.all.final)) < (0.2 * nrow(aucs.all.final))
    aucs.all.final <- aucs.all.final[, columns.to.keep]
    
    aucs.all.final[is.na(aucs.all.final)] <- 0
    
    aucs.all.final
}

ProjectFeatureDataPCA <- function(feature.data, num.components) {
    princ <- prcomp(feature.data)
    
    reduced.dim <- predict(princ, newdata=feature.data)[,1:num.components]
    colnames(feature.data) <- gsub("P", "F", colnames(feature.data))
    
    reduced.dim
}