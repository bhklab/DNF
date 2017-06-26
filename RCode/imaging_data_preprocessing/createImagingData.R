rm(list=ls())
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
#old_data <- readRDS("Data/combined_sens_datasets_with.RData")
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.names <- lincs.meta$pert_iname
lincs.names <- toupper(lincs.names)
lincs.names <- gsub(badchars, "", lincs.names)

imaging.meta <- read.delim("Data/Broad.HG005032.ProfilingData/imaging/cdrp.imaging.meta.cpd.txt", stringsAsFactors = FALSE)
imaging <- read.delim("Data/Broad.HG005032.ProfilingData/imaging/cdrp.imaging.profiles.txt", stringsAsFactors = FALSE)

# Clean names and IDs
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)
lincs.meta$pert_id <- toupper(lincs.meta$pert_id)
lincs.meta$pert_id <- gsub(badchars, "", lincs.meta$pert_id)
imaging.meta$name <- toupper(imaging.meta$name)
imaging.meta$name <- gsub(badchars, "", imaging.meta$name)

# Find relevant intersection of drugs in L1000 and cell painting. 
# This is based on matching the "name" column from the cell painting data
# against both the pert_iname and pert_id columns from the L1000 meta data,
# and then getting rid of dupes. 
pert.id.inter <- which(imaging.meta$name %in% lincs.meta$pert_id)
pert.name.inter <- which(imaging.meta$name %in% lincs.meta$pert_iname)
# Get rid of dupes
pert.name.inter <- imaging.meta$name[pert.name.inter]
pert.name.inter <- pert.name.inter[!duplicated(pert.name.inter)]
pert.name.inter <- match(pert.name.inter, imaging.meta$name)

combined <- unique(union(pert.id.inter, pert.name.inter))
# Get rid of dupes
temp.combined <- imaging.meta$name[combined]
temp.combined <- temp.combined[!duplicated(temp.combined)]
combined <- match(temp.combined, imaging.meta$name)

# This should be FALSE as we have gotten rid of dupes by this point
any(duplicated(imaging.meta[combined, "name"]))

imaging_drug_names <- imaging.meta[combined, c("BROAD_ID", "name")]
imaging_drug_names <- imaging_drug_names[!duplicated(imaging_drug_names$name),]
imaging.subsetted <- imaging[which(imaging$BROAD_ID %in% imaging_drug_names$BROAD_ID),]
imaging.subsetted <- imaging.subsetted[!duplicated(imaging.subsetted$BROAD_ID),]

aligned_indices <- match(imaging.subsetted$BROAD_ID, imaging_drug_names$BROAD_ID)
rownames(imaging.subsetted) <- imaging_drug_names$name[aligned_indices]
imaging.subsetted <- imaging.subsetted[, -1]

# Optional normalization
imaging.subsetted <- scale(imaging.subsetted)
imaging.subsetted <- imaging.subsetted[, colSums(is.na(imaging.subsetted)) != nrow(imaging.subsetted)]

# Optional dimensionality reduction
princ <- prcomp(imaging.subsetted, scale=TRUE)
n.comp <- 50

df.components <- predict(princ, newdata=imaging.subsetted)[,1:n.comp]

imaging.subsetted <- t(df.components)
imaging.subsetted <- as.matrix(imaging.subsetted)


saveRDS(imaging.subsetted, "Data/imaging_subsetted_predictive_pca.RData")
