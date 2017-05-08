library(PharmacoGx)

load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")
#load("PSets/CTRPv2.RData")
#load("PSets/GDSC1000.RData")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

datasets <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)
#datasets <- list(GDSC1000, CTRPv2)
cellines <- sort(unique(unlist(sapply(datasets, function(x) { return (cellNames(x)) }))))
drugs <- unique(unlist(sapply(datasets, function(x) { return (paste(pSetName(x), gsub(badchars, "",toupper(drugNames(x))), sep=":::")) })))

aucs.all <- matrix(NA, nrow=length(drugs), ncol=length(cellines), dimnames=list(drugs, cellines))

for (i in 1:length(datasets)) {
    aucs <- summarizeSensitivityProfiles(pSet = datasets[[i]], sensitivity.measure = "auc_recomputed")
    ### Clean up drug names the same way they're being cleaned in the DNF code. By doing this here
    ### we end up with drug rownames and colnames remaining in sorted order. Otherwise if we did this step
    ### in the DNF code, it might result in names that are out of order
    rownames(aucs) <- toupper(rownames(aucs))
    rownames(aucs) <- gsub(badchars, "", rownames(aucs))
    aucs.all[paste(pSetName(datasets[[i]]), rownames(aucs), sep=":::"), colnames(aucs)] <- aucs
}

### Calculate pearson correlation and zero out NA's. Note that zero 
### could be replaced with a more clever prior in the future.
aucs.cor <- cor(x=t(aucs.all), method="pearson", use="pairwise.complete.obs")
aucs.cor <- apply(aucs.cor, 1, function(x) ifelse(is.na(x),0,x))

### merge correlations for identical drugs across datasets
aucs.drugs <- sapply(strsplit(rownames(aucs.cor), ":::"), function (x) { return (x[[2]])})
aucs.drugs2 <- sort(unique(aucs.drugs))

### Create placeholder for final correlation matrix. Copy values from first correlation
### matrix for drugs that don't have duplicates
aucs.cor2 <- matrix(NA, nrow=length(aucs.drugs2), ncol=length(aucs.drugs2), dimnames=list(aucs.drugs2, aucs.drugs2))
iix <- sapply(strsplit(rownames(aucs.cor)[!duplicated(aucs.drugs)], ":::"), function (x) { return (x[[2]])})
aucs.cor2[iix, iix] <- aucs.cor[!duplicated(aucs.drugs), !duplicated(aucs.drugs)]

### Compute number of samples used for correlations in aucs.all. This is slightly trickier
### due to the fact that use="pairwise.complete.obs" is being used
num_samples_used = matrix(0, nrow=base::nrow(aucs.all), ncol=base::nrow(aucs.all))
dimnames(num_samples_used) <- dimnames(aucs.cor)
for (i in 1:base::nrow(aucs.all)) {
    for (j in 1:base::nrow(aucs.all)) {
        x = aucs.all[i,]
        len_x = length(x[!is.na(x)])
        
        y = aucs.all[j, ]
        len_y = length(y[!is.na(y)])
            
        num_samples_used[i, j] = min(len_x, len_y)
    }
}

### Compute fisher transformation of correlations in aucs.all
fisher_transformed <- fisherz(aucs.cor)

### Compute standard errors for aucs.all
standard_errors <- 1 / sqrt(num_samples_used - 3)
dimnames(standard_errors) <- dimnames(num_samples_used)


aucs.dupl <- aucs.drugs[duplicated(aucs.drugs)]
### Iterate over drugs that are found across multiple datasets, combine the corresponding
### Fisher-Z transformed correlations with survcomp::combine.est, and convert the combined value
### back into a correlation by using the inverse Fisher-Z transformation
self_discrepancies <- list()
for (i in 1:length(aucs.dupl)) {
    iix <- paste(sapply(datasets, pSetName), aucs.dupl[i], sep=":::")
    iix <- intersect(iix, rownames(aucs.all))
    self_discrepancies[[i]] <- aucs.cor[iix, iix]
    
    ### Currently not doing anything about a drug with itself. I.e. drug1 in dataset 1 and 2,
    ### only the correlation with drug 2 to drug n are being modified, but the correlation 
    ### between drug 1 and itself (which is normall 1 if it's only in 1 dataset), is being left
    ### alone. 
    indexed_rows <- fisher_transformed[iix, -which(colnames(fisher_transformed) %in% iix)]
    for (j in 1:ncol(indexed_rows)) {
        combined <- combine.est(fisher_transformed[iix, j], standard_errors[iix, j])
        drug_name <- strsplit(colnames(fisher_transformed)[j], ":::")[[1]][2]
        aucs.cor2[aucs.dupl[i], drug_name] <- fisherz(combined$estimate, inv=TRUE)
        aucs.cor2[drug_name, aucs.dupl[i]] <- fisherz(combined$estimate, inv=TRUE)
    }
}

dataset_pairs <- list()
### Iterate over self_discrepancies which is a list of matrices where each cell in the matrix
### contains a correlation for a given drug between different datasets. 
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
                
                dataset_pairs[[pair_name]][[dname]] <- discrepancy_matrix[j, k]
            }
        }
    }
}

unlisted_num_samples <- c()
for (i in 1:length(dataset_pairs)) {
    d_set_pair <- dataset_pairs[[i]]
    
    for (j in 1:length(d_set_pair)) {
        first_d_set <- strsplit(names(dataset_pairs)[i], '-')[[1]][1]
        second_d_set <- strsplit(names(dataset_pairs)[i], '-')[[1]][2]
        
        d_name <- names(d_set_pair)[j]
        print(paste(first_d_set, d_name, sep=':::'))
        print(paste(second_d_set, d_name, sep=':::'))
        
        n_samples <- num_samples_used[paste(first_d_set, d_name, sep=':::'), paste(second_d_set, d_name, sep=':::')]
        unlisted_num_samples <- c(unlisted_num_samples, n_samples)
    }
}


d <- data.frame(x = unlist(dataset_pairs), 
                grp = rep(names(dataset_pairs), times = sapply(dataset_pairs, length)),
                num_samples <- unlisted_num_samples)
ggplot(d,aes(x = grp, y = x)) + geom_boxplot()

boxplot(dataset_pairs, ylab="Correlation")
stripchart(dataset_pairs, vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, col = 'green')
title("Correlation of same Drug With Itself Across Dataset Pairs")

best_drugs <- lapply(dataset_pairs, function(x) {sort(x, decreasing = TRUE)[1:5]})
worst_drugs <- lapply(dataset_pairs, function(x) {sort(x, decreasing = FALSE)[1:5]})

saveRDS(aucs.cor2, './Data/combined_sens.RData')



