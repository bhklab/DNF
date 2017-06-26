ConstImagingLayerFlexible <- function(imaging.data) {
    # Creates an affinity matrix for the imaging layer. This is based on a 
    # pearson correlation matrix calculated given the imaging features in 
    # imaging.data.
    #
    # Args:
    #   imaging.data: A matrix where columns are drugs are rows are imaging features
    #                 for those drugs.
    #
    # Returns:
    #   An affinity matrix for the imaging layer.
    imaging.cor <- cor(imaging.data, method="pearson", use="pairwise.complete.obs")
    imaging.cor <- apply(imaging.cor, 1, function(x) ifelse(is.na(x),0,x))
    
    imaging.aff <- SNFtool::affinityMatrix(1 - imaging.cor, 20, 0.5)    
}


