PrintCIndices <- function(c.indices) {
    cat("c.indices values from each layer vs. the benchmark: \n")
    
    for (i in 1:length(c.indices)) {
        cat(names(c.indices)[i])
        cat(": ")
        cat(c.indices[[i]])
        cat("\n")
    }
}