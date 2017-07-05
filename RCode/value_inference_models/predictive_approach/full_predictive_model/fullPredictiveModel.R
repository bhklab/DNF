# This script predicts the features for the perturbation layer given features from 
# the 2 other standard layers.

rm(list=ls())
library(caret)
library(missForest)
library(Matrix)
library(PharmacoGx)
library(doParallel)
library(org.Hs.eg.db)

source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/flexible_layers/structureDataFlexible.R")
source("RCode/meta_analysis/meta_analysis_helpers.R")
source("RCode/averageAucs.R")
source("RCode/value_inference_models/predictive_approach/full_predictive_model/fullPredictiveModelHelpers.R")
source("RCode/flexible_layers/drugTargetBenchFlexible.R")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")

### Create the required perturbation data
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

pert.file.name <- "Data/L1000_compound_signatures.RData"
df.pert <- PerturbationDataFlexible(pert.file.name, lincs.meta)
df.pert <- t(df.pert)
df.pert <- df.pert[rowSums(is.na(df.pert)) < (0.05 * ncol(df.pert)), ]
colnames(df.pert) <- gsub("-", "_", colnames(df.pert))
###

### Create the required structure data. Note that the fingerprints are converted into
### a matrix where there is a 1 in a column if that bit in the fingerprint is on
fingerprints <- StructureDataFlexible(lincs.meta)

c1 <- makeCluster(8)
registerDoParallel(c1)

res <- foreach(i=1:length(fingerprints)) %dopar% {
    placeholder <- numeric(1024)
    fp <- fingerprints[[i]]
    bits <- fp@bits
    
    placeholder[bits] <- 1
    placeholder
}

stopCluster(c1)

df.fingerprints <- matrix(unlist(res), nrow=length(res), byrow=T)
rownames(df.fingerprints) <- names(fingerprints)

###

### Read the imaging features from file
imaging.data <- readRDS("Data/imaging_processed/imaging_subsetted_predictive_pca.RData")
imaging.data <- t(imaging.data)

datasets <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)

### Create sensitivity layer features
res <- CreateSensitivityFeatures(datasets, badchars, rownames(df.pert))
overlap <- res$overlap

# Figure out which drugs in the overlap of layers that we're trying to predict with
# have drug target info available
benchmark <- DrugTargetBenchFlexible(overlap, "temp_gmt", use.ctrpv2=TRUE, use.clue=FALSE,
                                     use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE)
benchmark.drugs <- rownames(benchmark)

# Create training and test sets for the model. The test set will consist only of
# drugs that have drug target info available.
set.seed(317)
test.indices <- sample(1:length(benchmark.drugs), size=25)
test.drugs <- benchmark.drugs[test.indices]
train.drugs <- setdiff(overlap, test.drugs)

sens.data.train <- res$sens.data.final[train.drugs, ]
pert.data.train <- df.pert[train.drugs, ]
strc.data.train <- df.fingerprints[train.drugs, ]

sens.data.test <- res$sens.data.final[test.drugs, ]
pert.data.test <- df.pert[test.drugs, ]
strc.data.test <- df.fingerprints[test.drugs, ]
#imaging.data.train <- imaging.data[train.drugs, ]

prediction.matrix.train <- cbind(pert.data.train, sens.data.train, strc.data.train)
input.layers <- list(sensitivity=sens.data.train, strucutre=strc.data.train)
output.layer <- pert.data.train

# Create the models used to predict feature values using just the training set data
models <- CreateModels(input.layers, output.layer, prediction.matrix.train, model.name="gaussprLinear")

# Predict feaure values on the test set data
prediction.matrix.test <- cbind(sens.data.test, strc.data.test)
predicted.perturbation.features <- MakePredictions(models, test.drugs, colnames(pert.data.test)[1:1],
                                                   prediction.matrix.test)
