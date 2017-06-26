ComputeCIndex <- function (x, y) {
    thresholds <- sort(x, decreasing = TRUE)
    sensitivities <- integer(length(thresholds) + 1)
    fpr <- integer(length(thresholds) + 1)
    
    sensitivities[1] <- 0
    fpr[1] <- 1
    auc <- 0
    
    for (i in 1:length(thresholds)) {
        threshold <- thresholds[i]   
        
        temp <- integer(length(x))
        temp[x < threshold] <- 0
        temp[x >= threshold] <- 1
        
        truth <- factor(y)
        predicted <- factor(temp)
        sens <- sensitivity(predicted, truth, positive=1)
        specif <- specificity(predicted, truth, negative = 0)    
        
        sensitivities[i] <- sens
        fpr[i] <- 1- specif
        
        if (i > 1) {
            auc <- auc + (sens * (fpr[i] - fpr[i - 1]))    
        } 
    }
    
    res <- list(tpr=sensitivities, fpr=fpr, auc=auc)
    
    return(res)
}

RandomGuessing <- function(y) {
    x <- runif(length(y))
    thresholds <- sort(x, decreasing = TRUE)
    
    sensitivities <- integer(length(thresholds) + 1)
    fpr <- integer(length(thresholds) + 1)
    
    sensitivities[1] <- 0
    fpr[1] <- 1
    auc <- 0
    
    for (i in 1:length(thresholds)) {
        threshold <- thresholds[i]   
        
        temp <- integer(length(y))
        #temp[x < threshold] <- 0
        #temp[x >= threshold] <- 1
        
        indices <- sample(x = length(y), size = 351, replace=FALSE)
        temp[indices] <- 1
        
        truth <- factor(y)
        predicted <- factor(temp)
        sens <- sensitivity(predicted, truth, positive=1)
        specif <- specificity(predicted, truth, negative = 0)    
        
        sensitivities[i] <- sens
        fpr[i] <- 1- specif
        
        if (i > 1) {
            auc <- auc + (sens * (fpr[i] - fpr[i - 1]))    
        } 
    }
    
    res <- list(tpr=sensitivities, fpr=fpr, auc=auc)
    
    return(res)
    
}

FindBadPerformers <- function(x, y) {
    bad.performers <- integer(length(x))
    false.positives <- integer(length(x))
    false.negatives <- integer(length(x))
    true.positives <- integer(length(x))
    true.negatives <- integer(length(x))
    
    thresholds <- sort(x, decreasing = TRUE)
    
    for (i in 1:length(thresholds)) {
        threshold <- thresholds[i]
        
        temp <- integer(length(x))
        temp[x < threshold] <- 0
        temp[x >= threshold] <- 1
        
        res <- integer(length(x))
        res[temp != y] <- 1
        
        false.pos <- integer(length(x))
        false.pos[temp == 1 & y == 0] <- 1
        
        false.neg <- integer(length(x))
        false.neg[temp == 0 & y == 1] <- 1
        
        true.pos <- integer(length(x))
        true.pos[temp == 1 & y == 1] <- 1
        
        true.neg <- integer(length(x))
        true.neg[temp == 0 & y == 0] <- 1
         
        bad.performers <- bad.performers + res
        false.positives <- false.positives + false.pos
        false.negatives <- false.negatives + false.neg
        true.positives <- true.positives + true.pos
        true.negatives <- true.negatives + true.neg
    }
    
    list(bad.performers=bad.performers, false.positives=false.positives, 
         false.negatives=false.negatives, true.positives=true.positives,
         true.negatives=true.negatives)
}