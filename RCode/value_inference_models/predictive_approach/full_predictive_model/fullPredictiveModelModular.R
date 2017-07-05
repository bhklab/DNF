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

datasets <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)

### Create the required perturbation data
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

pert.file.name <- "Data/L1000_compound_signatures.RData"
df.pert <- CreatePerturbationFeatures(pert.file.name, badchars)
###

### Create the required structure data. Note that the fingerprints are converted into
### a matrix where there is a 1 in a column if that bit in the fingerprint is on
df.fingerprints <- CreateStructureFeatures(lincs.meta)
###

### Read the imaging features from file
imaging.data <- readRDS("Data/imaging_processed/imaging_subsetted_predictive_pca.RData")
imaging.data <- t(imaging.data)

### Create sensitivity layer features
res <- CreateSensitivityFeatures(datasets, badchars, rownames(df.pert))
overlap <- res$overlap
rm(CCLE, GDSC1000, gCSI, FIMM, CTRPv2, datasets)
gc()
# Figure out which drugs in the overlap of layers that we're trying to predict with
# have drug target info available
benchmark <- DrugTargetBenchFlexible(overlap, "temp_gmt", use.ctrpv2=TRUE, use.clue=FALSE,
                                     use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE)
benchmark.drugs <- rownames(benchmark)

# Create training and test sets for the model. The test set will consist only of
# drugs that have drug target info available.

input.layer.names <- c("sensitivity")
output.layer.name <- c("perturbation")
model.name <- "gaussprLinear"
test.size <- 25

c1 <- makeCluster(8)
registerDoParallel(c1)

res <- Main(input.layer.names=input.layer.names, output.layer.name=output.layer.name,
     model.name=model.name, test.size=test.size, sens.data.final=res$sens.data.final,
     pert.data.final=df.pert, strc.data.final=df.fingerprints, imaging.data.final=NULLm,
     overlap=res$overlap, benchmark.drugs=benchmark.drugs)

Main <- function(input.layer.names, output.layer.name, model.name, test.size,
                 sens.data.final=NULL, pert.data.final=NULL, strc.data.final=NULL, imaging.data.final=NULL,
                 overlap, benchmark.drugs) {
    set.seed(317)
    
    all.layer.names <- union(input.layer.names, output.layer.name)
    
    test.indices <- sample(1:length(benchmark.drugs), size=25)
    test.drugs <- benchmark.drugs[test.indices]
    train.drugs <- setdiff(overlap, test.drugs)
    
    sens.data.train <- NULL
    pert.data.train <- NULL
    strc.data.train <- NULL
    imaging.data.train <- NULL
    
    sens.data.test <- NULL
    pert.data.test <- NULL
    strc.data.test <- NULL
    imaging.data.test <- NULL
    
    if ('sensitivity' %in% all.layer.names) {
        sens.data.train <- sens.data.final[train.drugs, ]
        sens.data.test <- sens.data.final[test.drugs, ]
    }
    
    if ('structure' %in% all.layer.names) {
        strc.data.train <- strc.data.final[train.drugs, ]
        strc.data.test <- strc.data.final[test.drugs, ]
    }
    
    if ('perturbation' %in% all.layer.names) {
        pert.data.train <- pert.data.final[train.drugs, ]
        pert.data.test <- pert.data.final[test.drugs, ]
    }
    
    if ('imaging' %in% all.layer.names) {
        imaging.data.train <- imaging.data.final[train.drugs, ]
        imaging.data.test <- imaging.data.final[test.drugs, ]
    }
    
    all.layers.train <- list(sensitivity=sens.data.train, structure=strc.data.train,
                             perturbation=pert.data.train, imaging=imaging.data.train)
    all.layers.test <- list(sensitivity=sens.data.test, structure=strc.data.test,
                            perturbation=pert.data.test, imaging=imaging.data.test)
    
    input.layers.train <- all.layers.train[input.layer.names]
    output.layer.train <- all.layers.train[[output.layer.name]]
    
    input.layers.test <- all.layers.test[input.layer.names]
    output.layer.test <- all.layers.test[[output.layer.name]]
    
    prediction.matrix.train <- do.call(cbind, all.layers.train[c(output.layer.name, input.layer.names)])
    
    # Create the models used to predict feature values using just the training set data
    models <- CreateModels(input.layers.train, output.layer.train, prediction.matrix.train, model.name=model.name)
    
    # Predict feaure values on the test set data
    prediction.matrix.test <- do.call(cbind, all.layers.test[c(output.layer.name, input.layer.names)])
    predicted.perturbation.features <- MakePredictions(models, test.drugs, colnames(output.layer.test),
                                                       prediction.matrix.test)
    
    res <- list(models=models, predicted.perturbation.features=predicted.perturbation.features)
}
