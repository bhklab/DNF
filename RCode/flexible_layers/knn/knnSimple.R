rm(list=ls())

source("RCode/flexible_layers/knn/drugTargetsKNN.R")
source("RCode/flexible_layers/compConcordIndxFlexible.R")
source("RCode/flexible_layers/generateROCPlotFlexible.R")
source("RCode/predPerf.R")

integrated <- readRDS("integrated.RData")
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
common.drugs <- colnames(integrated)
gmt.file.name <- "temp_gmt"
use.ctrpv2 <- TRUE
use.clue <- FALSE
use.chembl <- FALSE
use.dbank <- FALSE
use.dtc <- FALSE

data.bench <- DrugTargetsKNN(common.drugs, 
                             gmt.file.name, use.ctrpv2=use.ctrpv2,
                             use.clue=use.clue, use.chembl=use.chembl,
                             use.dbank=use.dbank, use.dtc=use.dtc) # 141 x 141 drug-drug adjacency matrix --> 141
data.bench[] <- lapply(data.bench, as.character)

num.targets <- length(unique(data.bench$TARGET_NAME))
drug.names <- sort(unique(data.bench$MOLECULE_NAME))
target.names <- sort(unique(data.bench$TARGET_NAME))

integrated <- integrated[drug.names, drug.names]

# Counts is a matrix of the shape NUM_DRUGS X NUM_TARGETS.
# Each drug will have counts correpsonding to each possible target.
# Most of these counts will be 0.
counts <- matrix(data=0, nrow=length(drug.names), ncol=num.targets)
rownames(counts) <- drug.names
colnames(counts) <- sort(unique(data.bench$TARGET_NAME))

# Neighbours is a a matirx of the shape NUM_DRUGS X K. 
# For example, a row looks like: Drug X     34 21 ... 10 4.
k <- 5
neighbours <- matrix(0, nrow=length(drug.names), ncol=k)
rownames(neighbours) <- drug.names

weights <- matrix(0, nrow=length(drug.names), ncol=k)
rownames(weights) <- drug.names

# For each drug, obtain the indices of the k closest neighbours
# according to the integrated similarity measure. For some reason,
# using existing KNN packages gave different results (I guess they
# thought that rows corresponds to points and columns correspond to 
# dimensions).
for (i in 1:nrow(neighbours)) {
    row <- integrated[i, ]
    res <- sapply(sort(row, index.return=TRUE), '[')
    res.indices <- res[, "ix"]
    res.weights <- res[, 'x']
    
    res.indices <- res.indices[(length(res.indices) - (k)):(length(res.indices) - 1)]
    res.weights <- res.weights[(length(res.weights) - (k)):(length(res.weights) - 1)]
    neighbours[i, ] <- res.indices
    weights[i, ] <- res.weights
}

# Convert the indices of nearest neigbhours into the actual drug
# names of those neighbours.
neighbours <- apply(neighbours, 1, function(x) {drug.names[x]})
neighbours <- t(neighbours)
rownames(neighbours) <- drug.names

# Now if we want to weight the targets based on the similarities between drugs,
# we have to determine some kind of threshold. For example, if the similarity
# between drug X and drug Y is 0.003, and the similarity between drug X and drug
# Z is 0.001, we should definitely trust the targets coming from drug Y more than 
# those coming from drug Z.

# Iterate over each drug in the neighbours matrix. For each neighbour
# of a given drug, obtain that neighbour's drug targets and use those 
# targets to increment the counts matrix in the apporpriate location.
for (i in 1:nrow(neighbours)) {
    drug <- rownames(neighbours)[i]
    names.of.neighbours <- neighbours[i, ]
    
    for (j in 1:length(names.of.neighbours)) {
        neighbour.name <- names.of.neighbours[j]
        
        relevant.targets <- data.bench[data.bench$MOLECULE_NAME == neighbour.name, "TARGET_NAME"]
        counts[drug, relevant.targets] <- counts[drug, relevant.targets] + (1 * weights[i, j])
        #counts[drug, relevant.targets] <- counts[drug, relevant.targets] + 1
    }
}

top.predictions <- sapply(1:nrow(counts), function(i, target.names, drug.names) {
    x <- counts[i, ]
    index.of.max <- which.max(x)
    res <- list()
    res[[drug.names[i]]] <- x[index.of.max]
    names(res[[drug.names[i]]]) <- target.names[index.of.max]
    res
}, target.names=colnames(counts), drug.names=rownames(counts))

num.correct.predictions <- 0
target.confusion.matrix <- matrix(0, nrow=num.targets, ncol=num.targets)
colnames(target.confusion.matrix) <- target.names
rownames(target.confusion.matrix) <- target.names

for (i in 1:length(top.predictions)) {
    drug.name <- names(top.predictions)[i]
    target.name <- names(top.predictions[[i]])

    possible.true.targets <- data.bench[data.bench$MOLECULE_NAME == drug.name, "TARGET_NAME"]
    
    print(paste("Drug Name", drug.name, sep=": "))
    print(paste("Target Name", target.name, sep=": "))
    print(paste("Possible Targets", possible.true.targets, sep=": "))
    
    if (target.name %in% possible.true.targets) {
        num.correct.predictions <- num.correct.predictions + 1
        
        diag(target.confusion.matrix)[possible.true.targets] <- diag(target.confusion.matrix)[possible.true.targets] + 1
    } else {
        target.confusion.matrix[possible.true.targets, target.name] <- target.confusion.matrix[possible.true.targets, target.name] + 1
    }
}

n = sum(target.confusion.matrix) # number of instances
nc = nrow(target.confusion.matrix) # number of classes
diag = diag(target.confusion.matrix) # number of correctly classified instances per class 
rowsums = apply(target.confusion.matrix, 1, sum) # number of instances per class
colsums = apply(target.confusion.matrix, 2, sum) # number of predictions per class
p = rowsums / n # distribution of instances over the actual classes
q = colsums / n # distribution of instances over the predicted classes

accuracy = sum(diag) / n 
accuracy 

precision = diag / colsums 
recall = diag / rowsums 
f1 = 2 * precision * recall / (precision + recall) 
data.frame(precision, recall, f1) 

macroPrecision = mean(precision, na.rm = TRUE)
macroRecall = mean(recall[recall != 0], na.rm = TRUE)
macroF1 = mean(f1, na.rm = TRUE)
data.frame(macroPrecision, macroRecall, macroF1)
