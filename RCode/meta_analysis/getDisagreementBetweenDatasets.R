GetDisagreementsBetweenDatasets <- function(self.discrepancies) {
    # Iterate over self.discrepancies which is a list of matrices where each cell 
    # in the matrix contains a correlation for a given drug between different datasets.
    #
    # Args:
    #   self.discrepancies: A list containing matrices. Each matrix is essentially
    #                       a sub-matrix of aucs.cor and has correlations of a drug
    #                       with itself across the various datasets that it is found in.
    #
    # Returns:
    #   A list of pairs of datasets. Each value in the list is a vector and every element
    #   of that vector is a correlation of drug with itself between the two datasets in
    #   the pair.
    dataset.pairs <- list()
    
    for (i in 1:length(self.discrepancies)) {
        discrepancy.matrix <- self.discrepancies[[i]]
        
        for (j in 1:nrow(discrepancy.matrix)) {
            for (k in (j +1):ncol(discrepancy.matrix)) {
                if (i !=  k && k <= ncol(discrepancy.matrix)) {
                    ### Obtain the name of the pair of datasets. For example, CCLE and GDSC1000
                    rname <- strsplit(rownames(discrepancy.matrix)[j], ":::")[[1]][1]
                    cname <- strsplit(colnames(discrepancy.matrix)[k], ":::")[[1]][1]
                    
                    ### We don't care about the correlations between a dataset and itself since those
                    ### are trivially 0, so we skip this loop iteration.
                    if (rname == cname) {
                        next
                    }
                    
                    ### Get the drug name
                    dname <- strsplit(rownames(discrepancy.matrix)[j], ":::")[[1]][2]
                    
                    pair.name <- paste(rname, cname, sep="-")
                    
                    ### Add an empty list to the list of dataset pairs if the current dataset pair doesn't
                    ### already exist.
                    if (is.null(dataset.pairs[pair.name])) {
                        dataset.pairs[[pair.name]] <- list()
                    }
                    
                    ### Find the correlation value between drug x in dataset A and drug x in dataset B
                    dataset.pairs[[pair.name]][[dname]] <- discrepancy.matrix[j, k]
                }
            }
        }
    }
    
    dataset.pairs
}