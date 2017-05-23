constImagingLayer <- function(imaging.data) {
    imaging.cor <- cor(imaging.data, method="pearson", use="pairwise.complete.obs")
    imaging.cor <- apply(imaging.cor, 1, function(x) ifelse(is.na(x),0,x))
    
    imaging.aff <- SNFtool::affinityMatrix(1 - imaging.cor, 20, 0.5)    
}


