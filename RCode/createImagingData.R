old_data <- readRDS("Data/combined_sens_datasets_with.RData")

imaging_meta <- read.delim("Data/Broad.HG005032.ProfilingData/imaging/cdrp.imaging.meta.cpd.txt", stringsAsFactors = FALSE)
imaging <- read.delim("Data/Broad.HG005032.ProfilingData/imaging/cdrp.imaging.profiles.txt", stringsAsFactors = FALSE)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

imaging_drugs <- imaging_meta$name
imaging_drugs <- toupper(imaging_drugs)
imaging_drugs <- gsub(badchars, "", imaging_drugs)

imaging_drug_names <- imaging_meta[imaging_drugs %in% colnames(old_data), c("BROAD_ID", "name")]
imaging_drug_names <- imaging_drug_names[!duplicated(imaging_drug_names$name),]
imaging_subsetted <- imaging[which(imaging$BROAD_ID %in% imaging_drug_names$BROAD_ID),]
imaging_subsetted <- imaging_subsetted[!duplicated(imaging_subsetted$BROAD_ID),]

aligned_indices <- match(imaging_subsetted$BROAD_ID, imaging_drug_names$BROAD_ID)
rownames(imaging_subsetted) <- imaging_drug_names$name[aligned_indices]
imaging_subsetted <- imaging_subsetted[, -1]
imaging_subsetted <- as.matrix(imaging_subsetted)
imaging_subsetted <- t(imaging_subsetted)

saveRDS(imaging_subsetted, "Data/imaging_subsetted.RData")