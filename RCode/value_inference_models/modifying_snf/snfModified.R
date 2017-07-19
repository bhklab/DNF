rm(list=ls())

set.seed(9833)

library(reshape2)
library(caret)
library(SNFtool)
library(org.Hs.eg.db)
source("RCode/flexible_layers/sensitivityDataFlexible.R")
source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/flexible_layers/structureDataFlexible.R")

source("RCode/flexible_layers/constSensitivityLayerFlexible.R")
source("RCode/flexible_layers/constPerturbationLayerFlexible.R")
source("RCode/flexible_layers/constStructureLayerFlexible.R")

source("RCode/cindexComp2.R")
source("RCode/flexible_layers/compConcordIndxFlexible.R")
source("RCode/flexible_layers/generateDrugPairsFlexible.R")
source("RCode/flexible_layers/drugTargetBenchFlexible.R")
source("RCode/flexible_layers/compDrugTargetBenchmarkFlexible.R")
source("RCode/flexible_layers/stackedLayerAnalysisHelpers.R")
source("RCode/flexible_layers/printCIndices.R")
source("RCode/flexible_layers/printPVals.R")
source("RCode/flexible_layers/generateROCPlotFlexible.R")
source("RCode/predPerf.R")
source("RCode/value_inference_models/modifying_snf/snfModifiedHelpers.R")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

pert.file.name <- "Data/L1000_compound_signatures.RData"
sensitivity.file.name <- "Data/combined_sensitivity//combined_sens_iname_replaced.RData"

sens.data <- SensitivityDataFlexible(sensitivity.file.name)  ## 645 X 239
pert.data <- PerturbationDataFlexible(pert.file.name, lincs.meta)  ## 978 X 23

layers <- list(sens.names = sort(colnames(sens.data)), pert.names=colnames(pert.data))
common.drugs <- Reduce(intersect, Filter(Negate(is.null),layers))
print(length(common.drugs))
sens.names <- rownames(sens.data)

lincs.meta.subset <- lincs.meta[match(sens.names, lincs.meta$pert_iname),]
lincs.meta.subset <- lincs.meta.subset[!is.na(lincs.meta.subset$X),]

strc.data <- StructureDataFlexible(lincs.meta.subset)  ## a vector  --> 239 elemnts
strc.cor <- fingerprint::fp.sim.matrix(strc.data, method = "tanimoto")
colnames(strc.cor) <- names(strc.data)
rownames(strc.cor) <- names(strc.data)

pert.cor <- cor(pert.data[, common.drugs], method = "pearson", use = "pairwise.complete.obs")

sens.cor <- sens.data

### Now make pert and strc layers as big as senitvity
missing.strc.drugs <- setdiff(sens.names, common.drugs)
missing.pert.drugs <- setdiff(sens.names, common.drugs)

pert.cor.augmented <- matrix(0, nrow=nrow(sens.cor), ncol=ncol(sens.cor))
rownames(pert.cor.augmented) <- rownames(sens.cor)
colnames(pert.cor.augmented) <- colnames(sens.cor)
pert.cor.augmented[rownames(pert.cor), colnames(pert.cor)] <- pert.cor

strc.cor.augmented <- matrix(0, nrow=nrow(sens.cor), ncol=ncol(sens.cor))
rownames(strc.cor.augmented) <- rownames(sens.cor)
colnames(strc.cor.augmented) <- colnames(sens.cor)
strc.cor.augmented[rownames(strc.cor), colnames(strc.cor)] <- strc.cor

strc.aff.mat <- SNFtool::affinityMatrix(1 - strc.cor.augmented, 20, 0.5)
sens.aff.mat <- SNFtool::affinityMatrix(1 - sens.cor, 20, 0.5)
pert.aff.mat <- SNFtool::affinityMatrix(1 - pert.cor.augmented, 20, 0.5)

strc.aff.mat[missing.strc.drugs, ] <- sens.aff.mat[missing.strc.drugs, ]
strc.aff.mat[, missing.strc.drugs] <- sens.aff.mat[, missing.strc.drugs]

pert.aff.mat[missing.pert.drugs, ] <- sens.aff.mat[missing.pert.drugs, ]
pert.aff.mat[, missing.pert.drugs] <- sens.aff.mat[, missing.pert.drugs]

integrated <- SNFtool::SNF(list(sens.aff=sens.aff.mat, strc.aff=strc.aff.mat, pert.aff=pert.aff.mat))
colnames(integrated) <- colnames(sens.aff.mat)
rownames(integrated) <- rownames(sens.aff.mat)

CompDrugTargetBenchmarkFlexible(common.drugs=common.drugs, gmt.file.name="temp",
                                strc.aff.mat=strc.aff.mat, sens.aff.mat=sens.aff.mat, pert.aff.mat=pert.aff.mat,
                                integration=integrated, luminex.aff.mat=NULL,
                                imaging.aff.mat=NULL, target.roc.file.name="blessed_roc.pdf",
                                use.ctrpv2=TRUE, use.clue=FALSE, use.chembl=FALSE, 
                                use.dbank=FALSE, use.dtc=FALSE)

### Now to come up with a general procedure for creating augmented matrices 
all.drugs <- sort(Reduce(union, list(rownames(strc.cor), rownames(sens.cor), rownames(pert.cor))))
correlation.matrices <- list(sens=sens.cor, pert=pert.cor, strc=strc.cor)
# augmented.matrices <- list(sens=NULL, pert=NULL, strc=NULL)

# 
# for (i in 1:length(augmented.matrices)) {
#     layer.name <- names(augmented.matrices)[i]
#     
#     augmented.matrices[[layer.name]] <- matrix(0, nrow=length(all.drugs), ncol=length(all.drugs))
#     colnames(augmented.matrices[[layer.name]]) <- all.drugs
#     rownames(augmented.matrices[[layer.name]]) <- all.drugs
# }
# 
# for (k in 1:length(augmented.matrices)) {
#     layer.name <- names(augmented.matrices)[k]
#     
#     drug.names <- rownames(correlation.matrices[[layer.name]])
#     
#     augmented.matrices[[layer.name]][drug.names, drug.names] <- correlation.matrices[[layer.name]]
# }
# 
# affinity.matrices <- list(sens=NULL, pert=NULL, strc=NULL)
# 
# for (i in 1:length(affinity.matrices)) {
#     layer.name <- names(affinity.matrices)[i]
#     
#     affinity.matrices[[layer.name]] <- SNFtool::affinityMatrix(1 - augmented.matrices[[layer.name]], 20, 0.5)
# }
# 
# for (i in 1:length(all.drugs)) {
#     drug.name <- all.drugs[i]
#     
#     for (k in 1:length(affinity.matrices)) {
#         layer.name <- names(affinity.matrices)[k]
#         other.layers <- names(affinity.matrices)[-k]
#         
#         drugs.not.in.layer <- setdiff(all.drugs, colnames(correlation.matrices[[layer.name]]))
#         
#         if (drug.name %in% rownames(correlation.matrices[[layer.name]])) {
#             for (d in drugs.not.in.layer) {
#                 values.from.other.layers <- vapply(other.layers, function(x) {
#                     if (d %in% rownames(correlation.matrices[[x]])) {
#                         return(affinity.matrices[[x]][drug.name, d])
#                     } else {
#                         return(NA)
#                     }
#                     
#                 }, numeric(1))                
#                 
#                 
#                 values.from.other.layers <- na.omit(values.from.other.layers)
#                 
#                 imputed.value <- mean(values.from.other.layers)
#                 
#                 affinity.matrices[[layer.name]][drug.name, d] <- imputed.value
#                 affinity.matrices[[layer.name]][d, drug.name] <- imputed.value
#             }            
#         } else {
#             for (d in all.drugs) {
#                 values.from.other.layers <- vapply(other.layers, function(x) {
#                     if (d %in% rownames(correlation.matrices[[x]])) {
#                         return(affinity.matrices[[x]][drug.name, d])
#                     } else {
#                         return(NA)
#                     }
#                     
#                 }, numeric(1))                
#                 
#                 
#                 values.from.other.layers <- na.omit(values.from.other.layers)
#                 
#                 imputed.value <- mean(values.from.other.layers)
#                 
#                 affinity.matrices[[layer.name]][drug.name, d] <- imputed.value
#                 affinity.matrices[[layer.name]][d, drug.name] <- imputed.value
#                 
#             }
#         }
#     }
# }

augmented.matrices <- CreateAugmentedMatrixSkeletons(all.drugs)
augmented.matrices <- ReplaceAugmentedExistingValues(augmented.matrices, correlation.matrices)
affinity.matrices <- CreateAffinityMatrices(augmented.matrices)
affinity.matrices <- ReplaceAffinityMatrixValues(affinity.matrices, correlation.matrices, all.drugs)

integrated <- SNFtool::SNF(affinity.matrices)
rownames(integrated) <- all.drugs
colnames(integrated) <- all.drugs

CompDrugTargetBenchmarkFlexible(common.drugs=common.drugs, gmt.file.name="temp",
                                strc.aff.mat=affinity.matrices$strc, 
                                sens.aff.mat=affinity.matrices$sens, pert.aff.mat=affinity.matrices$pert,
                                integration=integrated, luminex.aff.mat=NULL,
                                imaging.aff.mat=NULL, target.roc.file.name="blessed_roc.pdf",
                                use.ctrpv2=TRUE, use.clue=FALSE, use.chembl=FALSE, 
                                use.dbank=FALSE, use.dtc=FALSE)


