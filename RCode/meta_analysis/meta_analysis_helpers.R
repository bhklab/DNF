CreateAucsAll <- function(datasets, drugs, cell.lines, badchars) {
    # Creates a matrix whose rows are drugs and whose columns are cell lines. This is done by mergin
    # the different datasets such that drug names are prepended by the dataset that they came from, 
    # but cell lines aren't.
    #
    # Args:
    #   datasets: A list of PharmacoGx PSets.
    #   drugs: A vector of drug names prepended with the dataset they are from. The names are in 
    #          order with the dataset names in datasets.
    #   cell.lines: A vector of unique cell line names extracted from each PSet in datasets.
    #   badchars: A vector of characters to be removed from every drug name.
    #
    # Returns:
    #   A matrix whose rows are drugs and columns are cell lines. This can be used to 
    #   compute correlation between drugs.
    #
    aucs.all <- matrix(NA, nrow=length(drugs), ncol=length(cell.lines), dimnames=list(drugs, cell.lines))
    
    for (i in 1:length(datasets)) {
        aucs <- summarizeSensitivityProfiles(pSet = datasets[[i]], sensitivity.measure = "auc_recomputed")
        aucs <- aucs[!grepl(":", rownames(aucs)),]
        # Clean up drug names the same way they're being cleaned in the DNF code. By doing this here
        # we end up with drug rownames and colnames remaining in sorted order. Otherwise if we did this step
        # in the DNF code, it might result in names that are out of order
        old.names <- rownames(aucs)
        old.names <-  toupper(old.names)
        old.names <- gsub(badchars, "", old.names)
        
        # Use the drugs_with_ids.csv file to correct drug names. For example, abiraterone might be listed
        # using its IUPAC name instead in the PSet, so use that IUPAC name to find its common drug name.
        corrected.names <- ReplaceDrugNames(rownames(aucs), pSetName(datasets[[i]]))
        rownames(aucs) <- corrected.names
        rownames(aucs) <- toupper(rownames(aucs))
        rownames(aucs) <- gsub(badchars, "", rownames(aucs))
        
        # Replace with pert_iname from LINCS where possible
        rownames(aucs) <- ReplaceDrugNamesManualCurationPertInames(rownames(aucs))
        
        index.names <- paste(pSetName(datasets[[i]]), old.names, sep=":::")
        rownames(aucs.all)[which(rownames(aucs.all) %in% index.names)] <- paste(pSetName(datasets[[i]]), rownames(aucs), sep=":::")
        aucs.all[paste(pSetName(datasets[[i]]), rownames(aucs), sep=":::"), colnames(aucs)] <- aucs
    }
    
    aucs.all
}

ReplaceDrugNames <- function(drug.names, dataset) {
    # Uses the drugs_With_ids.csv file to replace the unique.id with the dataset.id.
    #
    # Args:
    #   drug.names: A vector of drug names without the PSet name prepended. These correspond to unique.ids
    #               that will be replaced by dataset.ids
    #   dataset: The name of the PSet in question.
    #
    # Returns:
    #   A vector of drug names that have been replaced based on what their names are in 
    #   the given dataset.
    if (dataset == "FIMM") {
        return(drug.names)
    }
    
    if (dataset == "gSCI") {
        dataset = "gCSI"
    }
    
    mapping.table <- read.csv("./Data/drugs_with_ids.csv", header=TRUE, stringsAsFactors = FALSE)
    corrected.names <- drug.names
    
    for (i in 1:length(drug.names)) {
        drug.name = drug.names[i]
        relevant.row <- mapping.table[which(mapping.table[,"unique.drugid"] == drug.name),]
        
        replace.value <- relevant.row[1, paste(dataset, "drugid", sep = ".")]
        
        if (!is.na(replace.value)) {
            corrected.names[i] <- replace.value        
        }
    }
    
    corrected.names
    
}

ReplaceDrugNamesManualCurationPertInames <- function(drug.names) {
    # Uses the file that Nehme manually curated to replace drug names
    # with their corresponding pert_inames from the L1000 dataset. The hope is
    # that this will improve overlap between the sensitivity and perturbation
    # layers.
    # 
    # Args:
    #   drug.names: A vector of drug names without the PSet name prepended. These
    #               will be used to index the drugnames column of the manually 
    #               curated file.
    #
    # Returns:
    #   A vector of drug names, some of which have been replaced by the 
    #   correponding pert_iname of the L1000 dataset.
    #
    manual.mapping <- read.csv("Data/mapping_manual_curation_inchikeys.csv")
    manual.mapping$drugnames <- toupper(manual.mapping$drugnames)
    manual.mapping$drugnames <- gsub(badchars, "", manual.mapping$drugnames)
    manual.mapping$pert_iname <- toupper(manual.mapping$pert_iname)
    manual.mapping$pert_iname <- gsub(badchars, "", manual.mapping$pert_iname)
    
    pert.inames <- drug.names
    
    for (i in 1:length(drug.names)) {
        drug.name <- drug.names[i]
        
        relevant.row <- manual.mapping[which(manual.mapping$drugnames == drug.name), ]
        
        replace.value <- relevant.row[1, "pert_iname"]
        
        if (!is.na(replace.value)) {
            pert.inames[i] <- replace.value
        }
    }
    
    pert.inames
}

CreateFinalCorrelationWithoutDups <- function(aucs.cor, aucs.drugs) {
    # Creates what will eventually be the final correlation matrix, except that
    # duplicated drugs have yet to have their correlations combined.
    #
    # Args:
    #   aucs.cor: A correlation matrix created using ComputeDatasetWideCorrelation where
    #             drug names are prepended with the dataset that they came from.
    #   aucs.drugs: The rownames of aucs.cor without the dataset preprended to them.
    #
    # Returns: 
    #   A correlation matrix without duplicated drug names. The correlations for 
    #   duplicated drugs have not yet been combined at this point.
    aucs.drugs2 <- sort(unique(aucs.drugs))
    
    ### Create placeholder for final correlation matrix. Copy values from first correlation
    ### matrix for drugs that don't have duplicates
    aucs.cor2 <- matrix(NA, nrow=length(aucs.drugs2), ncol=length(aucs.drugs2), dimnames=list(aucs.drugs2, aucs.drugs2))
    iix <- sapply(strsplit(rownames(aucs.cor)[!duplicated(aucs.drugs)], ":::"), function (x) { return (x[[2]])})
    aucs.cor2[iix, iix] <- aucs.cor[!duplicated(aucs.drugs), !duplicated(aucs.drugs)]
    
    aucs.cor2
}

ComputeNumSamplesUsed <- function(aucs.all, aucs.cor) {
    # Computes the number of samples used for each correlation in aucs.cor. The 
    # need for this comes from the fact that use="pairwise.complete.obs" is used
    # when computing aucs.cor.
    #
    # Args:
    #   aucs.all: A matrix created using the CreateAucsAll function where rownames are of the 
    #             form dataset_name:::drug_name and colnames are of the form cell_line_name.
    #   aucs.cor: The correlation matrix computed from aucs.all where both rownames 
    #             and column names take the form dataset_name:::drug_name.
    #
    # Returns:
    #   A matrix where each entry corresponds to the number of samples used to compute a
    #   correlation in aucs.cor at the same location. The rownames and colnames of this 
    #   matrix are identical to those of aucs.cor.
    num.samples.used = matrix(0, nrow=base::nrow(aucs.all), ncol=base::nrow(aucs.all))
    dimnames(num.samples.used) <- dimnames(aucs.cor)
    for (i in 1:base::nrow(aucs.all)) {
        for (j in 1:base::nrow(aucs.all)) {
            x = aucs.all[i,]
            x.len = length(x[!is.na(x)])
            
            y = aucs.all[j, ]
            y.len = length(y[!is.na(y)])
            
            num.samples.used[i, j] = min(x.len, y.len)
        }
    }
    
    num.samples.used
}

GetNumSamplesForDatasetPairs <- function(pair.name, dataset.pairs, num.samples.used) {
    # Indexes num.samples.used to obtain the number of samples used
    # in calculating the correlations for the drugs in dataset.pairs. 
    # This is used for coloring the correlation discrepancy boxplot based
    # on number of samples used.
    #
    # Args:
    #   pair.name: The concatenated result of two dataset names.
    #   dataset.pairs: 
    #   
    relevant.pair <- dataset.pairs[[pair.name]]
    result <- c()
    for (j in 1:length(relevant.pair)) {
        first.dataset <- strsplit(pair.name, '-')[[1]][1]
        second.dataset <- strsplit(pair.name, '-')[[1]][2]
        
        drug.name <- names(relevant.pair)[j]
        
        n.samples <- num.samples.used[paste(first.dataset, drug.name, sep=":::"), paste(second.dataset, drug.name, sep=":::")]
        result <- c(result, n.samples)
        names(result)[length(result)] <- paste(first.dataset, second.dataset, drug.name, sep="-")
    }
    
    result
}

CombineDrugCorrelations <- function(aucs.all, aucs.cor, aucs.cor2, aucs.dupl, datasets, fisher.transformed, standard.errors) {
    # Combines the correlations for duplicated drugs into a single value. For example
    # if drug x is found in 3 different datasets, and we want to know the correlation
    # it with drug y, we combine these 3 correlations into 1 using a Fisher Z transform.
    # 
    # Args:
    #   aucs.all: A matrix created by using CreateAucsAll(). Rownames are of the form dataset_name:::drug_name
    #             and colnames are of the form cell_line_name.
    #   aucs.cor: The original correlation matrix whose rownames and colnames are 
    #             of the form dataset_name:::drug_name. This will be used to determine the 
    #             so-called self discrepancies which means the same drug correlated with itself 
    #             in different datasets. If these correlations aren't high, this is a sign that 
    #             perhaps the methods used in the various assays were different.
    #   aucs.cor2: The final correlation matrix whose rownames and colnames are of the form drug_name. This is
    #              what will actually be modified and later returned. I.e. the values in fisher.transformed will
    #              be combined to obtain a new correlation.
    #   aucs.dupl: A vector of drugs that are found in mulitple datasets.
    #   datasets: A list of PharmacoGx PSets
    #   fisher.transformed: A matrix corresponding to the Fisher Z transformation applied to the correlations in
    #                       aucs.cor.
    #   standard.errors: A matrix containing the standard errors for all the Fisher Z statistics found 
    #                    fisher.transformed.
    #
    # Returns:
    #   A list of 2 values. One of the values is called self.discrepancies and is itself a list of 
    #   matrices where each matrix is a symmetric matrix containing correlations of the same drug with
    #   itself across various datasets. This is used later on to create a boxplot showing
    #   agreeement/disagreement between dataset pairs in terms of how correlated drug x
    #   in dataset 1 is with drug x in dataset 2. The second value, called aucs.cor2, is a modified
    #   version of aucs.cor2 where the correlations between drug x and drug y across various datasets 
    #   have been combined using combine.est.
    self.discrepancies <- list()
    for (i in 1:length(aucs.dupl)) {
        iix <- paste(sapply(datasets, pSetName), aucs.dupl[i], sep=":::")
        iix <- intersect(iix, rownames(aucs.all))
        self.discrepancies[[i]] <- aucs.cor[iix, iix]
        
        ### Currently not doing anything about a drug with itself. I.e. drug1 in dataset 1 and 2,
        ### only the correlation with drug 2 to drug n are being modified, but the correlation 
        ### between drug 1 and itself (which is normall 1 if it's only in 1 dataset), is being left
        ### alone. 
        indexed.rows <- fisher.transformed[iix, -which(colnames(fisher.transformed) %in% iix)]
        grouped.col.names <- list()
        
        splitted.names <- vapply(colnames(indexed.rows), strsplit, split=":::", FUN.VALUE = list(1))
        splitted.names <- do.call(rbind.data.frame, splitted.names)
        splitted.names[,1] <- as.character(splitted.names[,1])
        splitted.names[,2] <- as.character(splitted.names[,2])
        
        for (drug.name in unique(splitted.names[, 2])) {
            if (is.null(grouped.col.names[drug.name])) {
                grouped.col.names[[drug.name]] <- c()
            }
            
            grouped.col.names[[drug.name]] <- paste(splitted.names[which(splitted.names[, 2] == drug.name),1],
                                                    splitted.names[which(splitted.names[, 2] == drug.name),2],
                                                    sep=":::")
        }
        
        for (j in 1:length(grouped.col.names)) {
            col.names <- grouped.col.names[[j]]
            combined <- combine.est(as.vector(fisher.transformed[iix, col.names]), as.vector(standard.errors[iix, col.names]))
            drug.name <- strsplit(col.names[1], ":::")[[1]][2]
            aucs.cor2[aucs.dupl[i], drug.name] <- fisherz(combined$estimate, inv=TRUE)
            aucs.cor2[drug.name, aucs.dupl[i]] <- fisherz(combined$estimate, inv=TRUE)
        }
    }
    
    list(self.discrepancies = self.discrepancies, aucs.cor2 = aucs.cor2)
}

CombineDiscrepancyCorrelations <- function(aucs.all, aucs.cor, aucs.cor2, aucs.dupl, datasets, fisher.transformed, standard.errors) {
    # Combines the correlations of a drug with itself across multiple datasets. For example,
    # suppose that drug x is found in dataset 1, dataset 2, and dataset 3, the 
    # correlation is NOT trivially 1. Instead. the three different between-dataset correlations
    # need to be combined. This is done using the Fishze-Z transformation.
    #
    #
    # Args:
    #   aucs.all: A matrix created by using CreateAucsAll(). Rownames are of the form dataset_name:::drug_name
    #             and colnames are of the form cell_line_name.
    #   aucs.cor: The original correlation matrix whose rownames and colnames are 
    #             of the form dataset_name:::drug_name. This will be used to determine the 
    #             so-called self discrepancies which means the same drug correlated with itself 
    #             in different datasets. If these correlations aren't high, this is a sign that 
    #             perhaps the methods used in the various assays were different.
    #   aucs.cor2: The final correlation matrix whose rownames and colnames are of the form drug_name. This is
    #              what will actually be modified and later returned. I.e. the values in fisher.transformed will
    #              be combined to obtain a new correlation.
    #   aucs.dupl: A vector of drugs that are found in mulitple datasets.
    #   datasets: A list of PharmacoGx PSets
    #   fisher.transformed: A matrix corresponding to the Fisher Z transformation applied to the correlations in
    #                       aucs.cor.
    #   standard.errors: A matrix containing the standard errors for all the Fisher Z statistics found 
    #                    fisher.transformed.
    #
    # Returns:
    #   A matrix where all correlations have now been combined, not just the ones off the 
    #   diagonal. Note that this technically is no longer a correlation matrix since the 
    #   terms on the diagonal are no longer 1. This result has to be renormalized.
    #
    for (i in 1:length(aucs.dupl)) {
        iix <- paste(sapply(datasets, pSetName), aucs.dupl[i], sep=":::")
        iix <- intersect(iix, rownames(aucs.all))
        
        drug.name <- strsplit(iix[1], ":::")[[1]][2]
        
        upper.tri.indices <- upper.tri(fisher.transformed[iix, iix], diag = FALSE)
        relevant.fisher.vals <- fisher.transformed[iix, iix][upper.tri.indices]
        upper.tri.indices <- upper.tri(standard.errors[iix, iix], diag = FALSE)
        relevant.standard.errors <- standard.errors[iix, iix][upper.tri.indices]
        
        upper.tri.indices <- upper.tri(aucs.cor[iix, iix], diag = FALSE)
        relevant.correlations <- aucs.cor[iix, iix][upper.tri.indices]
        
        
        combined <- combine.est(relevant.fisher.vals, relevant.standard.errors)
        drug.name <- strsplit(colnames(fisher.transformed[iix, iix])[1], ":::")[[1]][2]
        aucs.cor2[drug.name, drug.name] <- fisherz(combined$estimate, inv=TRUE)
        
        #combined = weighted.mean(relevant.correlations, relevant.standard.errors)
        #aucs.cor2[drug.name, drug.name] <- combined
        
    
    }
    
    aucs.cor2
}

RenormalizeCorrelationMatrix <- function(aucs.cor2) {
    # Takes a pseudo-correlation matrix where the diagonal terms are not necessarily
    # 1 and does the following normalization: divide each row by the maximum term in that 
    # row, then simply paste 1's along the diagonal.
    # 
    # Args:
    #   aucs.cor2: A pseudo-correlation matrix (1 where the diagonal elements are
    #              not necessarily 1).
    #
    # Returns: A normalized version of the given correlation matrix.
    #
    for (i in 1:nrow(aucs.cor2)) {
        max.val <- max(aucs.cor2[i,])
        
        aucs.cor2[i,] = aucs.cor2[i,] / max.val
    }
    
    diag(aucs.cor2) <- 1
    
    aucs.cor2
}
