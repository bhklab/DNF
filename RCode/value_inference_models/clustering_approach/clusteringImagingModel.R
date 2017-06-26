rm(list=ls())
library(caret)
library(missForest)
library(Matrix)
library(PharmacoGx)
library(doParallel)
library(fields)
library(org.Hs.eg.db)

source("RCode/value_inference_models/clustering_approach/clusteringImagingModelHelpers.R")
source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/meta_analysis/meta_analysis_helpers.R")
source("RCode/averageAucs.R")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")

lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

pert.file.name <- "Data/L1000_compound_signatures.RData"
pert.data <- PerturbationDataFlexible(pert.file.name, lincs.meta)  ## 978 X 239
pert.data <- t(pert.data)

imaging.data <- readRDS("Data/imaging_subsetted_predictive_pca.RData")
imaging.data <- t(imaging.data)

datasets <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)

cell.lines <- sort(unique(unlist(sapply(datasets, function(x) { return (cellNames(x)) }))))
drugs <- unique(unlist(sapply(datasets, function(x) { return (paste(pSetName(x), gsub(badchars, "",toupper(drugNames(x))), sep=":::")) })))

### Merge datasets into one large matrix where rows are drugs (prepended with the corresponding dataset),
### and columns are cell lines
aucs.all <- CreateAucsAll(datasets, drugs, cell.lines = cell.lines, badchars = badchars)
rm(datasets)
rm(CCLE)
rm(GDSC1000)
rm(FIMM)
rm(gCSI)
rm(CTRPv2)
gc()

splitted <- strsplit(rownames(aucs.all), split=":::")
sens.names <- sapply(splitted, function(x) x[2])
overlap <- sort(intersect(rownames(imaging.data),
                          sens.names))

row.names.without.pset <- strsplit(rownames(aucs.all), ":::")
row.names.without.pset <- unique(sapply(row.names.without.pset, function(x){x[2]}))
# Create placeholder for final AUC matrix
aucs.all.final <- matrix(0, ncol=ncol(aucs.all), nrow=length(row.names.without.pset))
rownames(aucs.all.final) <- row.names.without.pset
colnames(aucs.all.final) <- colnames(aucs.all)

# Combine AUCs for drugs found in multiple datasets by averaging values. Also, only take subset of columns
# having less than 20% NA values
aucs.all.final <- AverageAUCS(row.names.without.pset, aucs.all, aucs.all.final)
aucs.all.final <- aucs.all.final[, colSums(is.na(aucs.all.final)) != nrow(aucs.all.final)]
columns.to.keep <- colSums(is.na(aucs.all.final)) < (0.2 * nrow(aucs.all.final))
aucs.all.final <- aucs.all.final[, columns.to.keep]

# Impute remaining missing value
registerDoParallel(6)
aucs.all.final <- missForest(aucs.all.final, maxiter=8, parallelize = "variables")$ximp

pert.data <- pert.data[intersect(rownames(pert.data), rownames(aucs.all.final)),]


names.of.drugs.to.augment <- sort(setdiff(sens.names, overlap))
drugs.to.augment <- matrix(0, ncol=ncol(imaging.data), nrow=length(names.of.drugs.to.augment))
rownames(drugs.to.augment) <- names.of.drugs.to.augment
colnames(drugs.to.augment) <- colnames(imaging.data)

drugs.to.augment <- PredictValuesBasedOnNeighbours(aucs.all.final, pert.data, imaging.data,
                                                   drugs.to.augment)

drugs.to.augment <- na.omit(drugs.to.augment)

final.result <- rbind(imaging.data, drugs.to.augment)
final.result <- t(final.result)

saveRDS(final.result, "imaging_data_cluster_pred.RData")
