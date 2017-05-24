PrintPVals <- function(p.vals) {
    cat("p-vals from the c.index comparison of integration layer vs. \n")
    
    for (i in 1:length(p.vals)) {
        cat(names(p.vals)[i])
        cat(": ")
        cat(p.vals[[i]])
        cat("\n")
    }
}