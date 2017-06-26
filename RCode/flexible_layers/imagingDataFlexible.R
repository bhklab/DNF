ImagingDataFlexible <- function(badchars) {
    # Loads in imaging data from file and cleans the drug names
    # by capitalizing them and removing bad chars
    #
    # Args:
    #   badchars: A string of bad characters to be removed from drug names.
    #
    # Returns:
    #   A matrix where each column corresponds to a drug, and each row
    #   is an imaging feature.
    imaging.data <- readRDS("Data/imaging_processed/imaging_subsetted_predictive_pca.RData")
    
    colnames(imaging.data) <- toupper(colnames(imaging.data))
    colnames(imaging.data) <- gsub(badchars, "", colnames(imaging.data))
    
    imaging.data <- imaging.data[, !duplicated(colnames(imaging.data))]
    imaging.data <- imaging.data[, order(colnames(imaging.data))]
    
    imaging.data
}