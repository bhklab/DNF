lincs <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
sens_data <- readRDS("Data/combined_sens_datasets_with.RData")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

lincs$pert_iname <- toupper(lincs$pert_iname)
lincs$pert_iname <- gsub(badchars, "", lincs$pert_iname)

intrsctLincsCombined <- intersect(lincs$pert_iname, rownames(sens_data))
lincsInters <- lincs[lincs$pert_iname %in% intrsctLincsCombined,,drop=F]
lincsInters <- lincsInters[!duplicated(lincsInters$pert_iname),,drop=F]
lincsInters$inchi_key

relevant_inchi_keys <- vapply(lincsInters$inchi_key, function(x) {
    strsplit(x, "=")[[1]][2]}, character(1))
result <- data.frame(inchi_keys=relevant_inchi_keys, SMILES=lincsInters$canonical_smiles)
View(result)
write.csv(result, "relevant_inchi_keys_and_smiles.csv", row.names=FALSE)
temp <- read.csv("relevant_inchi_keys_and_smiles.csv")
View(temp)
