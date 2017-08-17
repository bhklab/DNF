CalculateRankChanges <- function(neighbours.1, neighbours.2) {
    baseline <- 1:(ncol(neighbours.1))
    
    positive.changes <- numeric(length(rownames(neighbours.1)))
    names(positive.changes) <- rownames(neighbours.1)
    
    negative.changes <- numeric(length(rownames(neighbours.1)))
    names(negative.changes) <- rownames(neighbours.1)
    
    for (drug in rownames(neighbours.1)) {
        n.1 <- neighbours.1[drug, ]
        n.2 <- neighbours.2[drug, ]
        
        matched <- match(n.2, n.1)
        
        temp.diff <- baseline - matched
        
        positive.changes[drug] <- mean(temp.diff[temp.diff >=0])
        
        negative.changes[drug] <- mean(abs(temp.diff[temp.diff <= 0]))
    }
    
    mean.positive.change <- mean(positive.changes, na.rm=T)
    mean.negative.change <- mean(negative.changes, na.rm=T)
    
    return(list(positive.changes=positive.changes, 
                negative.changes=negative.changes,
                mean.positive.change=mean.positive.change,
                mean.negative.change=mean.negative.change))
}

CalculateRankChangesSpecific <- function(neighbours.1, neighbours.2, k) {
    positive.changes <- numeric(length(rownames(neighbours.1)))
    names(positive.changes) <- rownames(neighbours.1)
    
    negative.changes <- numeric(length(rownames(neighbours.1)))
    names(negative.changes) <- rownames(neighbours.1)
    
    for (drug in rownames(neighbours.1)) {
        n.1 <- neighbours.1[drug, ]
        n.2 <- neighbours.2[drug, 1:k]
        
        baseline <- 1:k
        matched <- match(n.2, n.1)
        
        temp.diff <- baseline - matched
        
        positive.changes[drug] <- mean(temp.diff[temp.diff >=0])
        
        negative.changes[drug] <- mean(abs(temp.diff[temp.diff <= 0]))
    }
    
    mean.positive.change <- mean(positive.changes, na.rm=T)
    mean.negative.change <- mean(negative.changes, na.rm=T)
    
    return(list(positive.changes=positive.changes,
                negative.changes=negative.changes,
                mean.positive.change=mean.positive.change,
                mean.negative.change=mean.negative.change))
}