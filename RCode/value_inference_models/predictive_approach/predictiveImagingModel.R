rm(list=ls())
library(caret)
library(missForest)
library(Matrix)
library(PharmacoGx)
library(doParallel)

source("RCode/predictiveImagingModelHelpers.R")
source("RCode/meta_analysis_helpers.R")
source("RCode/averageAucs.R")

load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")

imaging.data <- readRDS("Data/imaging_subsetted_pred_pca.RData")
imaging.data <- t(imaging.data)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

datasets <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)

cell.lines <- sort(unique(unlist(sapply(datasets, function(x) { return (cellNames(x)) }))))
drugs <- unique(unlist(sapply(datasets, function(x) { return (paste(pSetName(x), gsub(badchars, "",toupper(drugNames(x))), sep=":::")) })))

### Merge datasets into one large matrix where rows are drugs (prepended with the corresponding dataset),
### and columns are cell lines
aucs.all <- CreateAucsAll(datasets, drugs, cell.lines = cell.lines, badchars = badchars)

splitted <- strsplit(rownames(aucs.all), split=":::")
sens.names <- sapply(splitted, function(x) x[2])
overlap <- sort(intersect(rownames(imaging.data),
                          sens.names))
# Create placeholder for final AUC matrix
aucs.all.final <- matrix(0, ncol=ncol(aucs.all), nrow=length(overlap))
rownames(aucs.all.final) <- overlap
colnames(aucs.all.final) <- colnames(aucs.all)

# Combine AUCs for drugs found in multiple datasets by averaging values. Also, only take subset of columns
# having less than 20% NA values
aucs.all.final <- AverageAUCS(overlap, aucs.all, aucs.all.final)
columns.to.keep <- colSums(is.na(aucs.all.final)) < (0.2 * nrow(aucs.all.final))
aucs.all.final <- aucs.all.final[, columns.to.keep]

# Impute remaining missing value
aucs.all.final <- missForest(aucs.all.final, maxiter = 6)$ximp

imaging.data.final <- imaging.data[overlap, ]
prediction.matrix <- cbind(imaging.data.final, aucs.all.final)
models <- list()

# Build a regression model for each image feature and store it in the models list
for (i in 1:ncol(imaging.data.final)) {
    target.name <- colnames(imaging.data.final)[i]
    
    target.name.temp <- paste(target.name, " ~ ", sep="")
    
    feature.names <- colnames(aucs.all.final)
    feature.names <- sapply(feature.names, function(x) {paste('`', x, '`', sep="")})
    
    prediction.formula <- as.formula(paste(target.name.temp, paste(feature.names, sep="", collapse = " + "), sep=""))
    
    models[[target.name]] <- train(prediction.formula, data=prediction.matrix, method="gaussprLinear",
                                   na.action=na.pass)
}

# Take the sensitivity drugs that didn't intersect witht the imaging data,
# create a template matrix based on them.
names.of.drugs.to.augment <- sort(setdiff(sens.names, overlap))
drugs.to.augment <- matrix(0, ncol=ncol(aucs.all), nrow=length(names.of.drugs.to.augment))
rownames(drugs.to.augment) <- names.of.drugs.to.augment
colnames(drugs.to.augment) <- colnames(aucs.all)

# Average the AUCs for the drugs that are to be augmented. Again, 1 drug can 
# be found in multiple sensitivity datasets, so we take the average.
drugs.to.augment <- AverageAUCS(names.of.drugs.to.augment, aucs.all, drugs.to.augment)

# Keep only the columns that were used previously to build the models
drugs.to.augment <- drugs.to.augment[, columns.to.keep]

# Build parallel backend
registerDoParallel(4)

# Impute missing data
drugs.to.augment <- missForest(drugs.to.augment, maxiter=5, parallelize = "variables")$ximp

# Create a template to hold inferred imaging values for the sensitivity drugs that 
# were not found in the imaging data.
augmented.result <- matrix(0, nrow=nrow(drugs.to.augment), ncol=ncol(imaging.data.final))
augmented.result <- as.data.frame(augmented.result)
colnames(augmented.result) <- colnames(imaging.data.final)
rownames(augmented.result) <- rownames(drugs.to.augment)

# Infer columns one by one using the previously built models
for (i in 1:length(models)) {
    target.name <- names(models)[i]
    print(target.name)
    augmented.result[, target.name] <- predict(models[[target.name]], drugs.to.augment)
}

new.imaging.data <- rbind(imaging.data, augmented.result)
new.imaging.data <- t(new.imaging.data)

saveRDS(new.imaging.data, "Data/imaging_augmented.RData")

