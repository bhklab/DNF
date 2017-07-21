rm(list=ls())

set.seed(9833)

library(reshape2)
library(caret)
library(SNFtool)
library(org.Hs.eg.db)
source("RCode/flexible_layers/sensitivityDataFlexible.R")
source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/flexible_layers/structureDataFlexible.R")
source("RCode/flexible_layers/imagingDataFlexible.R")

source("RCode/flexible_layers/constImagingLayerFlexible.R")
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
imaging.file.name <- "Data/imaging_processed/imaging_subsetted_normalized_pca.RData"

sens.data <- SensitivityDataFlexible(sensitivity.file.name)  ## 645 X 239
pert.data <- PerturbationDataFlexible(pert.file.name, lincs.meta)  ## 978 X 23
imaging.data <- ImagingDataFlexible(badchars)

imaging.data <- imaging.data[, colnames(imaging.data) %in% colnames(sens.data)]

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
imaging.cor <- ConstImagingLayerFlexible(imaging.data)

sens.cor <- sens.data

### Now to come up with a general procedure for creating augmented matrices 
all.drugs <- sort(Reduce(union, list(rownames(strc.cor), rownames(sens.cor), rownames(pert.cor))))
correlation.matrices <- list(sens=sens.cor, pert=pert.cor, strc=strc.cor)

augmented.matrices <- CreateAugmentedMatrixSkeletons(names(correlation.matrices), all.drugs)
augmented.matrices <- ReplaceAugmentedExistingValues(augmented.matrices, correlation.matrices)
affinity.matrices <- CreateAffinityMatrices(augmented.matrices)
affinity.matrices <- ReplaceAffinityMatrixValuesFast(affinity.matrices, correlation.matrices, all.drugs)

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


