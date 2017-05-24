badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
#old_data <- readRDS("Data/combined_sens_datasets_with.RData")
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.names <- lincs.meta$pert_iname
lincs.names <- toupper(lincs.names)
lincs.names <- gsub(badchars, "", lincs.names)


luminex_meta <- read.delim("Data/Broad.HG005032.ProfilingData/luminex/cdrp.l1000.meta.cpd.txt", stringsAsFactors = FALSE)
luminex <- read.delim("Data/Broad.HG005032.ProfilingData/luminex/cdrp.l1000.profiles.txt", stringsAsFactors = FALSE)

luminex_drugs <- luminex_meta$name
luminex_drugs <- toupper(luminex_drugs)
luminex_drugs <- gsub(badchars, "", luminex_drugs)

luminex_drug_names <- luminex_meta[luminex_drugs %in% lincs.names, c("BROAD_ID", "name")]
luminex_drug_names <- luminex_drug_names[!duplicated(luminex_drug_names$name),]
luminex_subsetted <- luminex[which(luminex$BROAD_ID %in% luminex_drug_names$BROAD_ID),]
luminex_subsetted <- luminex_subsetted[!duplicated(luminex_subsetted$BROAD_ID),]

aligned_indices <- match(luminex_subsetted$BROAD_ID, luminex_drug_names$BROAD_ID)
rownames(luminex_subsetted) <- luminex_drug_names$name[aligned_indices]
luminex_subsetted <- luminex_subsetted[, -1]
luminex_subsetted <- as.matrix(luminex_subsetted)
luminex_subsetted <- t(luminex_subsetted)

saveRDS(luminex_subsetted, "Data/luminex_subsetted.RData")
