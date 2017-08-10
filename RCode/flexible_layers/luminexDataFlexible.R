LuminexDataFlexible <- function(badchars, file.path) {
    luminex.data <- readRDS(file.path)
    
    colnames(luminex.data) <- toupper(colnames(luminex.data))
    colnames(luminex.data) <- gsub(badchars, "", colnames(luminex.data))
    
    luminex.data <- luminex.data[, !duplicated(colnames(luminex.data))]
    luminex.data <- luminex.data[, order(colnames(luminex.data))]
    
    luminex.data
}