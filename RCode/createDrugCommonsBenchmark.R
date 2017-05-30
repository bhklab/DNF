sensitivity_data <- readRDS("Data/combined_sens_datasets_with.RData")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
#### See the number of common drug targets

dtc <- read.csv("Data/DTC_data.csv", stringsAsFactors = FALSE)
dtc_drug_names <- dtc$compound_name
dtc_drug_names <- toupper(dtc_drug_names)
dtc_drug_names <- gsub(badchars, "", dtc_drug_names)
intersect(colnames(sensitivity_data), dtc_drug_names)

#### Get the drug targets for common drugs
matches <- match(colnames(sensitivity_data), dtc_drug_names)
matches <- matches[!is.na(matches)]
dtc_subsetted <- dtc[matches,]

#### Only keep strong binders
dtc_subsetted <- dplyr::filter(dtc_subsetted, standard_type=="IC50")

dtc_subsetted <- dplyr::filter(dtc_subsetted,
                               (standard_units == "nM" & standard_value < 10000) | (standard_units == "microM" & standard_value < 10))

result <- dplyr::select(dtc_subsetted, compound_name, gene_names, assay_description)

drgTargets <- strsplit(result$gene_names, split = ",")
final <- data.frame(MOLECULE_NAME = rep(result$compound_name, sapply(drgTargets, length)), TARGET_NAME = unlist(drgTargets), stringsAsFactors = FALSE)

write.csv(final, file="Data/dtcTargets_filtered.csv", row.names=FALSE)