### Analyzes AUC by looking at false positives and false negatives. This is limited
### since pairs of drugs have to be taken into consideration. There is functionality
### implemented to suggest how well individual drugs perform, but this is flawed
### since we originally consider pairs of drugs and end up putting the blame
### on both drugs in the pair.
library(compiler)
library(ROCR)
source("RCode/temp/computeCIndex.R")

FindBadPerformersCompiled <- cmpfun(FindBadPerformers)

pairs <- readRDS("pairs.RData")

x <- pairs$integrPairs$obs.integr
y <- pairs$benchPairs$bench

res <- FindBadPerformersCompiled(x, y)
bad.performers <- res$bad.performers
false.positives <- res$false.positives
false.negatives <- res$false.negatives
true.positives <- res$true.positives
true.negatives <- res$true.negatives

to.analyze <- list(bad.performers=bad.performers)

for (i in 1:length(to.analyze)) {
    performance.metric <- to.analyze[[i]]
    names(performance.metric) <- paste(pairs$integrPairs$Var1, pairs$integrPairs$Var2, sep=":")
    quantile.filtered <- performance.metric[performance.metric > quantile(performance.metric, 0.001)]
    
    
    splitted <- strsplit(names(quantile.filtered), ":")
    concatenated.1 <- unlist(lapply(splitted, function(x) {x[1]}))
    concatenated.2 <- unlist(lapply(splitted, function(x) {x[2]}))    
    
    all.drug.names <- unique(union(concatenated.1, concatenated.2))
    individual.drugs <- integer(length(all.drug.names))
    names(individual.drugs) <- all.drug.names
    num.of.pairs <-  integer(length(all.drug.names))
    names(num.of.pairs) <- all.drug.names
    
    final <- lapply(names(quantile.filtered), function(x, bad.perf) {
        pair <- bad.perf[x]
        splitted <- strsplit(x, ":")
        drug.1 <- splitted[[1]][1]
        drug.2 <- splitted[[1]][2]
        individual.drugs[drug.1] <<- individual.drugs[drug.1] + pair
        individual.drugs[drug.2] <<- individual.drugs[drug.2] + pair
        
        num.of.pairs[drug.1] <<- num.of.pairs[drug.1] + 1
        num.of.pairs[drug.2] <<- num.of.pairs[drug.2] + 1
    }, bad.perf=quantile.filtered)
    
    print(sort(individual.drugs, decreasing = TRUE))
    print("***************************************")
    print("***************************************")
    print("***************************************")
}

individual.confidences <- individual.drugs / (length(x) * num.of.pairs)

names(bad.performers) <-  paste(pairs$integrPairs$Var1, pairs$integrPairs$Var2, sep=":")
filtered.bad.performers <- bad.performers[bad.performers > quantile(bad.performers, 0.95)]
sort(filtered.bad.performers, decreasing=TRUE)[1:100]
sort(bad.performers, decreasing=FALSE)[1:500]
plot(density(bad.performers))

confidences <- 1 - (bad.performers / length(x))
temp.1 <- sort(confidences[endsWith(names(confidences), "AFATINIB")])
confidence.adjusted.weights <- confidences * x
names(confidence.adjusted.weights) <- names(bad.performers)
temp.2 <- sort(confidence.adjusted.weights[endsWith(names(confidence.adjusted.weights), "AFATINIB")])