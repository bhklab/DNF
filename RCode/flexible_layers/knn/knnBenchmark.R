source("RCode/flexible_layers/knn/drugTargetsKNN.R")
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

# Counts is a matrix of the shape NUM_DRUGS X NUM_TARGETS.
# Each drug will have counts correpsonding to each possible target.
# Most of these counts will be 0.
counts <- matrix(data=0, nrow=length(common.drugs), ncol=num.targets)
rownames(counts) <- common.drugs
colnames(counts) <- sort(unique(data.bench$TARGET_NAME))

# Neighbours is a a matirx of the shape NUM_DRUGS X K. 
# For example, a row looks like: Drug X     34 21 ... 10 4.
k <- 5
neighbours <- matrix(0, nrow=length(common.drugs), ncol=k)
rownames(neighbours) <- common.drugs

# For each drug, obtain the indices of the k closest neighbours
# according to the integrated similarity measure. For some reason,
# using existing KNN packages gave different results (I guess they
# thought that rows corresponds to points and columns correspond to 
# dimensions).
for (i in 1:nrow(neighbours)) {
    row <- integrated[i, ]
    res <- sapply(sort(row, index.return=TRUE), '[')
    res <- res[, "ix"]
    res <- res[(length(res) - (k)):length(res)]
    print(length(res))
    neighbours[i, ] <- res[-length(res)]
}

# Convert the indices of nearest neigbhours into the actual drug
# names of those neighbours.
neighbours <- apply(neighbours, 1, function(x) {common.drugs[x]})
neighbours <- t(neighbours)
rownames(neighbours) <- common.drugs

# Iterate over each drug in the neighbours matrix. For each neighbour
# of a given drug, obtain that neighbour's drug targets and use those 
# targets to increment the counts matrix in the apporpriate location.
for (i in 1:nrow(neighbours)) {
    drug <- rownames(neighbours)[i]
    names.of.neighbours <- neighbours[i, ]
    
    for (j in 1:length(names.of.neighbours)) {
        neighbour.name <- names.of.neighbours[j]
        
        relevant.targets <- data.bench[data.bench$MOLECULE_NAME == neighbour.name, "TARGET_NAME"]
        counts[drug, relevant.targets] <- counts[drug, relevant.targets] + 1
    }
}



###### Code not yet adapted to KNN approach
print("Benchmark obtained")
pairs.target <- GenerateDrugPairsFlexible(data.bench, strc.aff=strc.aff.mat, sens.aff=sens.aff.mat,
                                          pert.aff=pert.aff.mat, integration=integration,
                                          luminex.aff=luminex.aff.mat, imaging.aff=imaging.aff.mat)

print("Drug pairs obtained")
res.target <- CompConcordIndxFlexible(pairs.target)

print("Concordance Index computed")
PrintCIndices(res.target$c.index.list)
PrintPVals(res.target$p.vals.list)
saveRDS(res.target$bad.performers, "bad_performers.RData")
saveRDS(pairs.target, "pairs.RData")

GenerateROCPlotFlexible(pairs.target, target.roc.file.name, nrow(data.bench))
    
