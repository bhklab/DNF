ImagingData <- function(itersc, badchars) {
    imaging.data <- readRDS("Data/imaging_subsetted.RData")
    
    colnames(imaging.data) <- toupper(colnames(imaging.data))
    colnames(imaging.data) <- gsub(badchars, "", colnames(imaging.data))
    
    imaging.data <- imaging.data[, match(itersc$pert_iname, colnames(imaging.data))]
    imaging.data <- imaging.data[, !duplicated(colnames(imaging.data))]
    imaging.data <- imaging.data[, order(colnames(imaging.data))]
    
    imaging.data
}