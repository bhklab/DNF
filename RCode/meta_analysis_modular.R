library(PharmacoGx)
library(survcomp)
library(ggplot2)

source('./RCode/getDisagreementBetweenDatasets.R')
source('./RCode/visualizeDatasetDiscrepancies.R')
source('./RCode/meta_analysis_helpers.R')
load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

datasets_with_gdsc <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)
datasets_without_gdsc <- list(CCLE, gCSI, CTRPv2, FIMM)

datasets <- list(datasets_without_gdsc = datasets_without_gdsc, datasets_with_gdsc = datasets_with_gdsc)
for (i in 1:length(datasets)) {
    temp <- strsplit(names(datasets)[i], '_')[[1]]
    save_dir <- paste('combined_sens_', temp[1], '_', temp[2], '.RData', sep="")
    discrepancies <- main(datasets[[i]], save_dir = save_dir)
    print(discrepancies)
}

main <- function(datasets, save_dir) {
    cellines <- sort(unique(unlist(sapply(datasets, function(x) { return (cellNames(x)) }))))
    drugs <- unique(unlist(sapply(datasets, function(x) { return (paste(pSetName(x), gsub(badchars, "",toupper(drugNames(x))), sep=":::")) })))
    
    ### Merge datasets into one large matrix where rows are drugs (prepended with the corresponding dataset),
    ### and columns are cell lines
    aucs.all <- create_aucs_all(datasets, drugs, cell_lines = cellines, badchars = badchars)
    
    ### Calculate pearson correlation and zero out NA's. Note that zero 
    ### could be replaced with a more clever prior in the future.
    aucs.cor <- cor(x=t(aucs.all), method="pearson", use="pairwise.complete.obs")
    percentage_of_na <- sum(is.na(aucs.cor)) / (nrow(aucs.cor) * ncol(aucs.cor))
    print(paste("Percentage of NAs in correlation:", percentage_of_na, sep=" "))
    most_nas_for_any_drug <- max(colSums(is.na(aucs.cor)))
    print(paste("Most NAs for any drug:", most_nas_for_any_drug, sep=" "))
    aucs.cor <- apply(aucs.cor, 1, function(x) ifelse(is.na(x),0,x))
    
    ### Create a new matrix based on aucs.cor, except that now the drug names will no longer
    ### be prepended by the dataset, and row and column names will be unique and in sorted order
    ### Note that this matrix will have NAs for drugs that occur in multiple datasets, since the 
    ### correlations for these drugs must be merged later on using the Fisher Z stat method.
    ### *** The above sentence is not true. The way that the original aucs.cor is being indexed
    ### returns a matrix the full size of aucs.cor2. The role the duplicated() function plays is the following:
    ### for drugs that exist in multiple datasets, take the correlation that's in the first dataset.
    aucs.drugs <- sapply(strsplit(rownames(aucs.cor), ":::"), function (x) { return (x[[2]])})
    aucs.cor2 <- create_final_correlation_unmerged_dups(aucs.cor, aucs.drugs)
    
    ### Compute number of samples used for correlations in aucs.all. This is slightly trickier
    ### due to the fact that use="pairwise.complete.obs" is being used
    num_samples_used = compute_num_samples_used(aucs.all, aucs.cor)
    
    ### Compute fisher transformation of correlations in aucs.all
    fisher_transformed <- fisherz(aucs.cor)
    
    ### Compute standard errors for aucs.all
    standard_errors <- 1 / sqrt(num_samples_used - 3)
    dimnames(standard_errors) <- dimnames(num_samples_used)
    
    aucs.dupl <- aucs.drugs[duplicated(aucs.drugs)]
    ### Iterate over drugs that are found across multiple datasets, combine the corresponding
    ### Fisher-Z transformed correlations with survcomp::combine.est, and convert the combined value
    ### back into a correlation by using the inverse Fisher-Z transformation
    result <- combine_drug_correlations(aucs.all, aucs.cor, aucs.cor2, aucs.dupl, datasets, fisher_transformed, standard_errors)
    aucs.cor2 <- result$aucs.cor2
    self_discrepancies <- result$self_discrepancies
    
    ### Combine the drugs on the diagonal (i.e. same drug with itself across multiple datasets)
    # aucs.cor2 <- combine_disrepancy_correlations(aucs.all, aucs.cor2, aucs.dupl, datasets, fisher_transformed, standard_errors)
    
    dataset_pairs <- get_disagreements_between_datasets(self_discrepancies = self_discrepancies)
    unlisted_num_samples <- unlist(sapply(names(sapply(dataset_pairs, names)), get_num_samples_for_dataset_pairs, dataset_pairs = dataset_pairs, num_samples_used = num_samples_used))
    visualize_dataset_discrepancies(dataset_pairs = dataset_pairs, unlisted_num_samples = unlisted_num_samples)   
    
    saveRDS(aucs.cor2, save_dir)
    
    self_discrepancies
}


best_drugs <- lapply(dataset_pairs, function(x) {sort(x, decreasing = TRUE)[1:5]})
worst_drugs <- lapply(dataset_pairs, function(x) {sort(x, decreasing = FALSE)[1:5]})

saveRDS(aucs.cor2, './Data/combined_sens.RData')



