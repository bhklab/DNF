rm(list=ls())

library(MatrixCorrelation)
library(foreach)
library(doSNOW)
.libPaths('/mnt/work1/users/bhklab/Rlib/')
library(PharmacoGx)

c1 <- makeCluster(24, outfile="")
registerDoSNOW(c1)

signatures <- readRDS("signatures.RData")

concatenated.signatures <- list()

for (cell.line in names(signatures)) {
    cell.line.signatures <- signatures[[cell.line]][, , "estimate"]
    
    for (drug in colnames(cell.line.signatures)) {
        if (is.null(concatenated.signatures[[drug]])) {
            concatenated.signatures[[drug]] <- cell.line.signatures[, drug]
        } else {
            concatenated.signatures[[drug]] <- rbind(concatenated.signatures[[drug]], cell.line.signatures[, drug])
        }
    }
}

concatenated.signatures <- lapply(concatenated.signatures, 
function(x) {
  if (class(x) != "matrix") {
      x = as.matrix(x)
      return(t(x))
  } else {
      return(x)
  }
})
  

progress <- function(n) cat(sprintf("task %d is complete\n", n))
opts <- list(progress=progress)
  
correlation.matrix <- foreach(i=1:length(concatenated.signatures), 
                              .packages=c("MatrixCorrelation", "foreach")) %dopar% { 
    foreach (j=i:length(concatenated.signatures), .combine="c", .packages=c("MatrixCorrelation", "foreach")) %do% {
        require(foreach)
        d1 <- names(concatenated.signatures)[i]
        d2 <- names(concatenated.signatures)[j]
        
        mat.1 <- concatenated.signatures[[d1]]
        mat.2 <- concatenated.signatures[[d2]]
        temp <- RV(t(mat.1), t(mat.2))
        names(temp) <- d2
        temp
    }
}

names(correlation.matrix) <- names(concatenated.signatures)

saveRDS(correlation.matrix, "pert_matrix_correlations.rds")