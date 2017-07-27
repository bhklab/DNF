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

# correlation.matrix <- matrix(0, nrow=length(concatenated.signatures), ncol=length(concatenated.signatures))
# 
# rownames(correlation.matrix) <- names(concatenated.signatures)
# colnames(correlation.matrix) <- names(concatenated.signatures)

# for (i in names(concatenated.signatures)) {
#     mat.1 <- concatenated.signatures[[i]]
#     
#     for (j in names(concatenated.signatures)) {
#         if (j > i) {
#             mat.2 <- concatenated.signatures[[j]]
#             temp <- RV(t(mat.1), t(mat.2))
#             
#             correlation.matrix[i, j] <- temp
#             correlation.matrix[j, i] <- temp   
#         }
#     }
#     print(i)
# }

correlation.matrix <- foreach(i=names(concatenated.signatures), .combine="rbind") %dopar% 
{
    require(foreach)
    require(MatrixCorrelation)
    mat.1 <- concatenated.signatures[[i]]
    
    foreach (j=names(concatenated.signatures), .combine="c") %do% {
        mat.2 <- concatenated.signatures[[j]]
        temp <- RV(t(mat.1), t(mat.2))
        names(temp) <- j
        temp
    }
    
    print(i)
}

rownames(correlation.matrix) <- names(concatenated.signatures)

# diag(correlation.matrix) <- 1

saveRDS(correlation.matrix, "pert_matrix_correlations.rds")