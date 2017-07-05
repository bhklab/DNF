rm(list=ls())

source("RCode/flexible_layers/perturbationDataFlexible.R")
source("RCode/value_inference_models/predictive_approach/full_predictive_model/fullPredictiveModelHelpers.R")

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

### Create the required perturbation data
impute.lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
impute.lincs.meta$pert_iname <- toupper(impute.lincs.meta$pert_iname)
impute.lincs.meta$pert_iname <- gsub(badchars, "", impute.lincs.meta$pert_iname)

pert.file.name <- "Data/L1000_compound_signatures.RData"
impute.df.pert <- CreatePerturbationFeatures(pert.file.name, badchars)

predicted.perturbation.features <- readRDS("pred_pert_features_sens.RData")

df.pert.imputed <- impute.df.pert
df.pert.imputed[rownames(predicted.perturbation.features), ] <- predicted.perturbation.features
df.pert.imputed <- t(df.pert.imputed)

saveRDS(df.pert.imputed, "pert_data_imputed_sens.RData")
