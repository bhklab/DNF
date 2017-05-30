ConstLuminexLayerFlexible <- function(luminex.data) {
    luminex.cor <- cor(luminex.data, method="pearson", use="pairwise.complete.obs")
    luminex.cor <- apply(luminex.cor, 1, function(x) ifelse(is.na(x),0,x))
    
    luminex.aff <- SNFtool::affinityMatrix(1 - luminex.cor, 20, 0.5)    
}


