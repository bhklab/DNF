### SNF Explained

library(SNFtool)

features <- matrix(rnorm(5000), nrow=100, ncol=50)
similarities <- cor(features, method = "pearson")

# Based on analysis done below, it seems that affinity matrix represents the full kernel, perhaps not
# fully normalized.
aff.mat <- affinityMatrix(similarities)

rowSums(aff.mat)
View(aff.mat)

Dist1 <- dist2(as.matrix(Data1), as.matrix(Data1))
Dist2 <- dist2(as.matrix(Data2), as.matrix(Data2))
W1 <- affinityMatrix(Dist1)
W2 <- affinityMatrix(Dist2)

normalized.w1 <- W1 / rowSums(W1)
normalized.w1 <- (normalized.w1 + t(normalized.w1)) / 2

W1[, 101:200] <- 0
W1[101:200, ] <- 0

p.1 <- rbind(c(0.5, 0.2, 0.3), c(0.2, 0.6, 0.2), c(0.1, 0.2, 0.7))
p.2 <- rbind(c(0.8, 0.2, 0), c(0.3, 0.7, 0), c(0, 0, 0))

Wall <- list(W1, W2)

res.imputed <- tempSNF(Wall, 20, 20)

### .dominateSet source code
tempDominateSet <- function(xx,KK=20) {
    ###This function outputs the top KK neighbors.	
    
    zero <- function(x) {
        s = sort(x, index.return=TRUE)
        x[s$ix[1:(length(x)-KK)]] = 0
        return(x)
    }
    normalize <- function(X) X / rowSums(X)
    A = matrix(0,nrow(xx),ncol(xx));
    for(i in 1:nrow(xx)){
        A[i,] = zero(xx[i,]);
        
    }
    print(A)
    for (i in 1:nrow(A)) {
        if (sum(A[i, ]) == 0) {
            A[i, i] <- 1
        }
    }
    
    
    return(normalize(A))
}

### SNF source Code
tempSNF <- function (Wall, K = 20, t = 20) 
{
    LW = length(Wall)
    normalize <- function(X) X/rowSums(X)
    newW <- vector("list", LW)
    nextW <- vector("list", LW)
    for (i in 1:LW) {
        # The normalization being done here doesn't seem to correspond to what's happening 
        # in equation 2 of DNF. This is mainly due to the fact that row sums are not 1
        Wall[[i]] = normalize(Wall[[i]])
        Wall[[i]] = (Wall[[i]] + t(Wall[[i]]))/2
        temp <- Wall[[i]] 
        temp[is.na(temp)] <- 0
        Wall[[i]] <- temp
    }
    
    
    # As noted further below, newW is the sparse kernel. This suggests that the 
    # .dominateset is what further normalized Wall s.t. the row sums are all 1
    for (i in 1:LW) {
        newW[[i]] = (tempDominateSet(Wall[[i]], K))
    }
    for (i in 1:t) {
        for (j in 1:LW) {
            # Kind of redundant to index at [1] and [2] since the matrices in Wall should
            # all be square
            sumWJ = matrix(0, dim(Wall[[j]])[1], dim(Wall[[j]])[2])
            for (k in 1:LW) {
                if (k != j) {
                    sumWJ = sumWJ + Wall[[k]]
                }
            }
            # P^(v) = S^(v)(\frac{\sum\limits_{k \neq v} P^(k)}{m - 1})(S^(v))^T
            # This is equation 7 in the methods section of the SNF paper
            # So newW is the sparse kernel.
            nextW[[j]] = newW[[j]] %*% (sumWJ/(LW - 1)) %*% t(newW[[j]])
        }
        for (j in 1:LW) {
            Wall[[j]] = nextW[[j]] + diag(nrow(Wall[[j]]))
            Wall[[j]] = (Wall[[j]] + t(Wall[[j]]))/2
        }
    }
    W = matrix(0, nrow(Wall[[1]]), ncol(Wall[[1]]))
    for (i in 1:LW) {
        W = W + Wall[[i]]
    }
    W = W/LW
    W = normalize(W)
    W = (W + t(W) + diag(nrow(W)))/2
    return(W)
}
