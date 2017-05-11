### Iterate over self_discrepancies which is a list of matrices where each cell in the matrix
### contains a correlation for a given drug between different datasets.
get_disagreements_between_datasets <- function(self_discrepancies) {
    dataset_pairs <- list()
    
    for (i in 1:length(self_discrepancies)) {
        discrepancy_matrix <- self_discrepancies[[i]]
        
        for (j in 1:nrow(discrepancy_matrix)) {
            for (k in (j +1):ncol(discrepancy_matrix)) {
                if (i !=  k && k <= ncol(discrepancy_matrix)) {
                    ### Obtain the name of the pair of datasets. For example, CCLE and GDSC1000
                    rname <- strsplit(rownames(discrepancy_matrix)[j], ":::")[[1]][1]
                    cname <- strsplit(colnames(discrepancy_matrix)[k], ":::")[[1]][1]
                    
                    ### We don't care about the correlations between a dataset and itself since those
                    ### are trivially 0, so we skip this loop iteration.
                    if (rname == cname) {
                        next
                    }
                    
                    ### Get the drug name
                    dname <- strsplit(rownames(discrepancy_matrix)[j], ":::")[[1]][2]
                    
                    pair_name <- paste(rname, cname, sep="-")
                    
                    ### Add an empty list to the list of dataset pairs if the current dataset pair doesn't
                    ### already exist.
                    if (is.null(dataset_pairs[pair_name])) {
                        dataset_pairs[[pair_name]] <- list()
                    }
                    
                    ### Find the correlation value between drug x in dataset A and drug x in dataset B
                    dataset_pairs[[pair_name]][[dname]] <- discrepancy_matrix[j, k]
                }
            }
        }
    }
    
    dataset_pairs
}