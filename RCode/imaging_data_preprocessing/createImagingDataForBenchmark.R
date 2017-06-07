rm(list=ls())
source("RCode/flexible_layers/drugTargetBenchFlexible.R")
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
#old_data <- readRDS("Data/combined_sens_datasets_with.RData")

imaging.meta <- read.delim("Data/Broad.HG005032.ProfilingData/imaging/cdrp.imaging.meta.cpd.txt", stringsAsFactors = FALSE)
imaging <- read.delim("Data/Broad.HG005032.ProfilingData/imaging/cdrp.imaging.profiles.txt", stringsAsFactors = FALSE)

# Clean names and IDs
imaging.meta$name <- toupper(imaging.meta$name)
imaging.meta$name <- gsub(badchars, "", imaging.meta$name)

benchmark <- DrugTargetBenchFlexible(imaging.meta$name, "", use.ctrpv2=TRUE, use.clue=TRUE,
                                     use.chembl=TRUE, use.dbank=TRUE, use.dtc=TRUE)

imaging.meta <- imaging.meta[imaging.meta$name %in% colnames(benchmark), ]

imaging_drug_names <- imaging.meta[, c("BROAD_ID", "name")]
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
princ <- prcomp(imaging.subsetted)
n.comp <- 20

df.components <- predict(princ, newdata=imaging.subsetted)[,1:n.comp]
set.seed(111)
df.components <- df.components[sample(1:nrow(df.components), 350), ] # take a smaller sample since computationally expensive to do 1000

imaging.subsetted <- t(df.components)
imaging.subsetted <- as.matrix(imaging.subsetted)

saveRDS(imaging.subsetted, "Data/imaging_benchmark_pca.RData")
