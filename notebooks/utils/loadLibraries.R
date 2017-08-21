library(PharmacoGx)
library(apcluster)
library(rcdk)
library(fingerprint)
library(annotate)
library(org.Hs.eg.db)
library(SNFtool)
library(ROCR)
library(survcomp)
library(reshape2)
library(proxy)
library(PRROC)
library(apcluster)
library(profvis)

library(doParallel)
library(foreach)

source("RCode/flexible_layers/structureDataFlexible.R")
source("RCode/flexible_layers/sensitivityDataFlexible.R")
source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/flexible_layers/perturbationDataFlexibleCustomSig.R")
source("RCode/flexible_layers/luminexDataFlexible.R")
source("RCode/flexible_layers/imagingDataFlexible.R")

source("RCode/flexible_layers/constSensitivityLayerFlexible.R")
source("RCode/flexible_layers/constStructureLayerFlexible.R")
source("RCode/flexible_layers/constLuminexLayerFlexible.R")
source("RCode/flexible_layers/constImagingLayerFlexible.R")
source("RCode/flexible_layers/constPerturbationLayerFlexible.R")

source("RCode/flexible_layers/integrateLayersFlexible.R")
source("RCode/flexible_layers/drugTargetBenchFlexible.R")
source("RCode/flexible_layers/stackedLayerAnalysisHelpers.R")

source("RCode/cindexComp2.R")
source("RCode/predPerf.R")

source("RCode/knn/knnCV.R")
source("RCode/knn/knnHelpers.R")
source("RCode/knn/drugTargetsKNN.R")

source("RCode/goldenberg_imputation/medianSimilarity.R")
source("RCode/value_inference_models/modifying_snf/snfModifiedHelpers.R")

source("RCode/foreach_utils/foreachUtils.R")

# Load lincs metadata and clean the pert_inames in it
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"