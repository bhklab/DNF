LuminexDataFlexible <- function(badchars) {
    luminex.data <- readRDS("Data/luminex_subsetted.RData")
    
    colnames(luminex.data) <- toupper(colnames(luminex.data))
    colnames(luminex.data) <- gsub(badchars, "", colnames(luminex.data))
    
    luminex.data <- luminex.data[, !duplicated(colnames(luminex.data))]
    luminex.data <- luminex.data[, order(colnames(luminex.data))]
    
    luminex.data
}