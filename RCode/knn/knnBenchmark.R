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

integrated <- integrated[drug.names, drug.names]

# Counts is a matrix of the shape NUM_DRUGS X NUM_TARGETS.
# Each drug will have counts correpsonding to each possible target.
# Most of these counts will be 0.
counts <- matrix(data=0, nrow=length(drug.names), ncol=num.targets)
rownames(counts) <- drug.names
colnames(counts) <- sort(unique(data.bench$TARGET_NAME))

# Neighbours is a a matirx of the shape NUM_DRUGS X K. 
# For example, a row looks like: Drug X     34 21 ... 10 4.
k <- 9
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
    }
}

predicted.pairs <- melt(counts)
colnames(predicted.pairs)[1] <- "Var1"
colnames(predicted.pairs)[2] <- "Var2"
colnames(predicted.pairs)[3] <- "obs.integr"
predicted.indices <- which(predicted.pairs$obs.integr > 0)

bench.pairs <- matrix(0, nrow=length(drug.names), ncol=num.targets)
rownames(bench.pairs) <- drug.names
colnames(bench.pairs) <- sort(unique(data.bench$TARGET_NAME))

for (i in 1:nrow(bench.pairs)) {
    drug.name <- rownames(bench.pairs)[i]
    
    relevant.targets <- data.bench[data.bench$MOLECULE_NAME == drug.name, "TARGET_NAME"]
    
    bench.pairs[drug.name, relevant.targets] <- 1
}

bench.pairs <- melt(bench.pairs)
colnames(bench.pairs)[1] <- "Var1"
colnames(bench.pairs)[2] <- "Var2"
colnames(bench.pairs)[3] <- "bench"
bench.indices <- which(bench.pairs$bench > 0)

combined.indices <- union(predicted.indices, bench.indices)
predicted.pairs <- predicted.pairs[combined.indices, ]
bench.pairs <- bench.pairs[combined.indices, ]

pairs.list <- list()

pairs.list[["integrPairs"]] <- predicted.pairs
pairs.list[["benchPairs"]] <- bench.pairs

###### Code not yet adapted to KNN approach
print("Benchmark obtained")

print("Drug pairs obtained")
res.target <- CompConcordIndxFlexible(pairs.list)

perf.results <- predPerf(predicted.pairs$obs.integr, bench.pairs$bench)
f1.scores <- perf.results$f1.score@y.values[[1]]
f1.scores[1] <- 0
best.index <- which(f1.scores == max(f1.scores))

threshold <- perf.results$f1.score@x.values[[1]][best.index]

predictions <- predicted.pairs
predictions$obs.integr[predictions$obs.integr < threshold] <- 0
predictions$obs.integr[predictions$obs.integr >= threshold] <- 1


GenerateROCPlotFlexible(pairs.list, "temp_roc.pdf", length(drug.names))

