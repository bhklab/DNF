rm(list=ls())

source("RCode/flexible_layers/drugTargetBenchFlexible.R")
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
k <- 7
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
    
    max.x <- max(weights[i, ])
    fudge.factor <- - ((log(3)) / ((weights[i, ncol(weights) - 1] - weights[i, ncol(weights)]))) + 0.001
    sigmoid <- function(z) {1 / (1 + exp(- fudge.factor * (z - max.x)))}
    
    for (j in 1:length(names.of.neighbours)) {
        neighbour.name <- names.of.neighbours[j]
        
        relevant.targets <- data.bench[data.bench$MOLECULE_NAME == neighbour.name, "TARGET_NAME"]
        counts[drug, relevant.targets] <- counts[drug, relevant.targets] + (sigmoid(weights[i, j]))
        #counts[drug, relevant.targets] <- counts[drug, relevant.targets] + 1
    }
}

counts <- apply(counts, 1, function(x) {
    indices <- which(x > 0)
    res <- numeric(length(x))
    res[indices] <- exp(x[indices]) / sum(exp(x[indices]))
    res
})

counts <- t(counts)
colnames(counts) <- sort(unique(data.bench$TARGET_NAME))

top.predictions <- sapply(1:nrow(counts), function(i, target.names, drug.names) {
    x <- counts[i, ]
    sorted <- sort(x, index.return=T)
    index.of.predictions <- tail(sorted$ix[sorted$x > 0], 5)
    res <- list()
    res[[drug.names[i]]] <- x[index.of.predictions]
    names(res[[drug.names[i]]]) <- target.names[index.of.predictions]
    res
}, target.names=colnames(counts), drug.names=rownames(counts))

num.correct.predictions <- integer(nrow(counts))
names(num.correct.predictions) <- rownames(counts)
false.positive.targets <- integer(ncol(counts))
names(false.positive.targets) <- colnames(counts)

for (i in 1:length(top.predictions)) {
    drug.name <- names(top.predictions)[i]
    target.names <- names(top.predictions[[i]])
    
    possible.true.targets <- data.bench[data.bench$MOLECULE_NAME == drug.name, "TARGET_NAME"]
    
    print(paste("Drug Name", drug.name, sep=": "))
    print(paste("Predicted Targets", target.names, sep=": "))
    print(paste("Possible Targets", possible.true.targets, sep=": "))
    
    if (length(intersect(target.names, possible.true.targets)) > 0) {
        num.correct.predictions[drug.name] <- 1
    }
    
    false.positives <- setdiff(target.names, possible.true.targets)
    if (length(false.positives) > 0) {
        false.positive.targets[false.positives] <- false.positive.targets[false.positives] + 1
    }
}

accuracy <- sum(num.correct.predictions) / length(num.correct.predictions)

