create_aucs_all <- function(datasets, drugs, cell_lines, badchars) {
    aucs.all <- matrix(NA, nrow=length(drugs), ncol=length(cell_lines), dimnames=list(drugs, cell_lines))
    
    for (i in 1:length(datasets)) {
        aucs <- summarizeSensitivityProfiles(pSet = datasets[[i]], sensitivity.measure = "auc_recomputed")
        ### Clean up drug names the same way they're being cleaned in the DNF code. By doing this here
        ### we end up with drug rownames and colnames remaining in sorted order. Otherwise if we did this step
        ### in the DNF code, it might result in names that are out of order
        old_names <- rownames(aucs)
        old_names <-  toupper(old_names)
        old_names <- gsub(badchars, "", old_names)
        corrected_names <- replace_drug_names(rownames(aucs), pSetName(datasets[[i]]))
        rownames(aucs) <- corrected_names
        rownames(aucs) <- toupper(rownames(aucs))
        rownames(aucs) <- gsub(badchars, "", rownames(aucs))
        
        index_names <- paste(pSetName(datasets[[i]]), old_names, sep=":::")
        rownames(aucs.all)[which(rownames(aucs.all) %in% index_names)] <- paste(pSetName(datasets[[i]]), rownames(aucs), sep=":::")
        aucs.all[paste(pSetName(datasets[[i]]), rownames(aucs), sep=":::"), colnames(aucs)] <- aucs
    }
    
    aucs.all
}

replace_drug_names <- function(drug_names, dataset) {
    if (dataset == "FIMM") {
        return(drug_names)
    }
    
    if (dataset == "gSCI") {
        dataset = "gCSI"
    }
    
    mappingTable <- read.csv("./Data/drugs_with_ids.csv", header=TRUE, stringsAsFactors = FALSE)
    corrected_names <- drug_names
    
    for (i in 1:length(drug_names)) {
        d_name = drug_names[i]
        relevant_row <- mappingTable[which(mappingTable[,"unique.drugid"] == d_name),]
        
        replace_value <- relevant_row[1, paste(dataset, "drugid", sep = ".")]
        
        if (!is.na(replace_value)) {
            corrected_names[i] <- replace_value        
        }
    }
    
    corrected_names
    
}

### 
compute_correlation <- function(aucs) {
    aucs.cor <- cor(x=t(aucs.all), method="pearson", use="pairwise.complete.obs")
    aucs.cor <- apply(aucs.cor, 1, function(x) ifelse(is.na(x),0,x))
    
    aucs.cor
}

create_final_correlation_unmerged_dups <- function(aucs.cor, aucs.drugs) {
    ### 
    aucs.drugs2 <- sort(unique(aucs.drugs))
    
    ### Create placeholder for final correlation matrix. Copy values from first correlation
    ### matrix for drugs that don't have duplicates
    aucs.cor2 <- matrix(NA, nrow=length(aucs.drugs2), ncol=length(aucs.drugs2), dimnames=list(aucs.drugs2, aucs.drugs2))
    iix <- sapply(strsplit(rownames(aucs.cor)[!duplicated(aucs.drugs)], ":::"), function (x) { return (x[[2]])})
    aucs.cor2[iix, iix] <- aucs.cor[!duplicated(aucs.drugs), !duplicated(aucs.drugs)]
    
    aucs.cor2
}

compute_num_samples_used <- function(aucs.all, aucs.cor) {
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
    
    num_samples_used
}

get_num_samples_for_dataset_pairs <- function(pair_name, dataset_pairs, num_samples_used) {
    relevant_pair <- dataset_pairs[[pair_name]]
    result <- c()
    for (j in 1:length(relevant_pair)) {
        first_d_set <- strsplit(pair_name, '-')[[1]][1]
        second_d_set <- strsplit(pair_name, '-')[[1]][2]
        
        drug_name <- names(relevant_pair)[j]
        
        n_samples <- num_samples_used[paste(first_d_set, drug_name, sep=":::"), paste(second_d_set, drug_name, sep=":::")]
        result <- c(result, n_samples)
        names(result)[length(result)] <- paste(first_d_set, second_d_set, drug_name, sep="-")
    }
    
    result
}

combine_drug_correlations <- function(aucs.all, aucs.cor, aucs.cor2, aucs.dupl, datasets, fisher_transformed, standard_errors) {
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
        grouped_col_names <- list()
        
        splitted_names <- vapply(colnames(indexed_rows), strsplit, split=":::", FUN.VALUE = list(1))
        splitted_names <- do.call(rbind.data.frame, splitted_names)
        splitted_names[,1] <- as.character(splitted_names[,1])
        splitted_names[,2] <- as.character(splitted_names[,2])
        
        for (drug_name in unique(splitted_names[, 2])) {
            if (is.null(grouped_col_names[drug_name])) {
                grouped_col_names[[drug_name]] <- c()
            }
            
            grouped_col_names[[drug_name]] <- paste(splitted_names[which(splitted_names[, 2] == drug_name),1],
                                                    splitted_names[which(splitted_names[, 2] == drug_name),2],
                                                    sep=":::")
        }
        
        for (j in 1:length(grouped_col_names)) {
            col_names <- grouped_col_names[[j]]
            combined <- combine.est(as.vector(fisher_transformed[iix, col_names]), as.vector(standard_errors[iix, col_names]))
            drug_name <- strsplit(col_names[1], ":::")[[1]][2]
            aucs.cor2[aucs.dupl[i], drug_name] <- fisherz(combined$estimate, inv=TRUE)
            aucs.cor2[drug_name, aucs.dupl[i]] <- fisherz(combined$estimate, inv=TRUE)
        }
    }
    
    print(colnames(fisher_transformed))
    list(self_discrepancies = self_discrepancies, aucs.cor2 = aucs.cor2)
}

combine_disrepancy_correlations <- function(aucs.all, aucs.cor2, aucs.dupl, datasets, fisher_transformed, standard_errors) {
    for (i in 1:length(aucs.dupl)) {
        iix <- paste(sapply(datasets, pSetName), aucs.dupl[i], sep=":::")
        iix <- intersect(iix, rownames(aucs.all))
        
        upper_tri_indices <- upper.tri(fisher_transformed[iix, iix], diag = FALSE)
        relevant_fisher_Vals <- fisher_transformed[iix, iix][upper_tri_indices]
        upper_tri_indices <- upper.tri(standard_errors[iix, iix], diag = FALSE)
        relevant_standard_errors <- standard_errors[iix, iix][upper_tri_indices]
        
        combined <- combine.est(relevant_fisher_Vals, relevant_standard_errors)
        drug_name <- strsplit(colnames(fisher_transformed[iix, iix])[1], ":::")[[1]][2]
        
        aucs.cor2[drug_name, drug_name] <- fisherz(combined$estimate, inv=TRUE)
    
    }
    
    aucs.cor2
}

