old_data <- readRDS("Data/combined_sens_datasets_with.RData")

luminex_meta <- read.delim("Data/Broad.HG005032.ProfilingData/luminex/cdrp.l1000.meta.cpd.txt", stringsAsFactors = FALSE)
luminex <- read.delim("Data/Broad.HG005032.ProfilingData/luminex/cdrp.l1000.profiles.txt", stringsAsFactors = FALSE)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

luminex_drugs <- luminex_meta$name
luminex_drugs <- toupper(luminex_drugs)
luminex_drugs <- gsub(badchars, "", luminex_drugs)
luminex_drugs <- luminex_drugs[!duplicated(luminex_drugs)]

luminex_drug_names <- luminex_meta[which(luminex_drugs %in% colnames(old_data)), c("BROAD_ID", "name")]
luminex_subsetted <- luminex[which(luminex$BROAD_ID %in% luminex_drug_names$BROAD_ID),]
luminex_subsetted <- luminex_subsetted[!duplicated(luminex_subsetted$BROAD_ID),]

aligned_indices <- match(luminex_subsetted$BROAD_ID, luminex_drug_names$BROAD_ID)
rownames(luminex_subsetted) <- luminex_drug_names$name[aligned_indices]
luminex_subsetted <- luminex_subsetted[, -1]
luminex_subsetted <- as.matrix(luminex_subsetted)
luminex_subsetted <- t(luminex_subsetted)

saveRDS(luminex_subsetted, "Data/luminex_subsetted.RData")