rm(list=ls())
library(caret)
library(mice)
library(Matrix)
library(PharmacoGx)

source("RCode/meta_analysis_helpers.R")
source("RCode/averageAucs.R")

load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")

imaging.data <- readRDS("Data/imaging_subsetted_pca.RData")
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
# Index aucs.all using only drugs in overlap. Will be dupes since there are multiple datasets per drug.
# Remove dupes 
#aucs.all.final <- aucs.all[sens.names %in% overlap, ]

aucs.all.final <- matrix(0, ncol=ncol(aucs.all), nrow=length(overlap))
rownames(aucs.all.final) <- overlap
colnames(aucs.all.final) <- colnames(aucs.all)

aucs.all.final <- AverageAUCS(overlap, aucs.all, aucs.all.final)
aucs.all.final[is.na(aucs.all.final)] <- 0

#sens.names <- sens.names[sens.names %in% overlap]
#aucs.all.final <- aucs.all.final[!duplicated(sens.names),]

# Change the rownames of aucs.all to no longer have ::: in them
#splitted <- strsplit(rownames(aucs.all.final), split=":::")
#sens.names <- sort(sapply(splitted, function(x) x[2]))
#aucs.all.final <- aucs.all.final[sens.names, ]
#rownames(aucs.all.final) <- sens.names

imaging.data.final <- imaging.data[overlap, ]

prediction.matrix <- cbind(imaging.data.final, aucs.all.final)

models <- list()

for (i in 1:ncol(imaging.data.final)) {
    target.name <- colnames(imaging.data.final)[i]
    
    target.name <- paste(target.name, " ~ ", sep="")
    
    feature.names <- colnames(aucs.all.final)
    feature.names <- sapply(features.names, function(x) {paste('`', x, '`', sep="")})
    
    prediction.formula <- as.formula(paste(target.name, paste(feature.names, sep="", collapse = " + "), sep=""))
    
    models[[target.name]] <- train(prediction.formula, data=prediction.matrix, method="lm",
                                   na.action=na.pass)
}

new.drug <- aucs.all.final[1, ,drop=FALSE]
new.drug[is.na(new.drug)] <- 0
predict(models[[1]], new.drug)

