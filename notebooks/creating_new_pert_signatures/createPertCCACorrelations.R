rm(list=ls())

library(PMA)
library(foreach)
library(doParallel)
.libPaths('/mnt/work1/users/bhklab/Rlib/')
library(PharmacoGx)

c1 <- makeCluster(8, outfile="")
registerDoParallel(c1)

signatures <- readRDS("Data/signatures.RData")

concatenated.signatures <- list()

for (cell.line in names(signatures)) {
    cell.line.signatures <- signatures[[cell.line]][, , "estimate"]
    
    for (drug in colnames(cell.line.signatures)) {
        if (is.null(concatenated.signatures[[drug]])) {
            concatenated.signatures[[drug]] <- matrix(cell.line.signatures[, drug], nrow=1)
            colnames(concatenated.signatures[[drug]]) <- names(cell.line.signatures[, drug])
        } else {
            concatenated.signatures[[drug]] <- rbind(concatenated.signatures[[drug]], cell.line.signatures[, drug])
        }
        
        rownames(concatenated.signatures[[drug]])[nrow(concatenated.signatures[[drug]])] <- cell.line
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

correlation.matrix <- foreach(i=1:length(concatenated.signatures), 
                              .packages=c("CCA", "foreach")) %dopar% { 
                                  foreach (j=i:length(concatenated.signatures), .combine="c", .packages=c("CCA", "foreach")) %do% {
                                      require(foreach)
                                      d1 <- names(concatenated.signatures)[i]
                                      d2 <- names(concatenated.signatures)[j]
                                      
                                      mat.1 <- t(concatenated.signatures[[d1]])
                                      mat.2 <- t(concatenated.signatures[[d2]])
                                      
                                      intersection <- intersect(colnames(mat.1), colnames(mat.2))
                                      
                                      if (length(intersection) == 0) {
                                          temp <- NA
                                      } else {
                                          mat.1 <- mat.1[, intersection, drop=F]
                                          mat.2 <- mat.2[, intersection, drop=F]
                                          
                                          temp <- CCA::cc(mat.1, mat.2)$cor[1]
                                      }
                                      
                                      names(temp) <- d2
                                      temp
                                  }
                              }

names(correlation.matrix) <- names(concatenated.signatures)

saveRDS(correlation.matrix, "pert_matrix_cca_correlations.rds")
