#' Takes two numerical vectors and computes the concordance index between them
#' by comparing the order of values for two pairs of data each time
#'
#' This function return the concordance index and its p-value
#' along with the lower and upper confidence intervals of said p-value.
#'
#'
#' @examples
#' gdsc.auc <- summarizeSensitivityProfiles(GDSC, sensitivity.measure='auc_published')
#' xx <- summarizeSensitivityProfiles(predicted.by.model, gdsc.auc$Erlotinib)
#'
#' @param predictions {numeric} A vector of predicted drug responces which could be either continuous or discrete
#' @param observations {numeric} A vector of observed continuous drug responces
#' @param cutoff {numeric} A drug responce threshold which is used to classify cells to sensitive vs resistant to drug.
#' @param delta {numeric} The minimunm reliable difference between two drug sensitivity values to be considered as significantly various responses.
#' default value for delta is picked by looking into delta auc values between biological replicates across three
#' large pharmacogenomic studies, CTRPv2(370 drugs over ~15-20 cells) , GDSC(1 drug over ~600 cells), GRAY (85 drugs over ~10-50)
#' @param alpha {numeric} alpha level to compute confidence interval
#' @param outx {boolean} set to TRUE to not count pairs of observations tied on x as a relevant pair.
#' This results in a Goodman-Kruskal gamma type rank correlation.
#' @param alternative {character} what is the alternative hypothesis? Must be one of "two.sides", "less", and "greater".
#' @return [list] ! list of concordance index and its pvalue
#' along with the lower and upper confidence intervals
#' @export

PairedConcordanceIndex <- function(predictions, observations, delta.pred=0.2, delta.obs=0.2, alpha = 0.05, outx = TRUE, alternative = c("two.sided", "less", "greater"),  logic.operator=c("or", "and")) {
    alternative <- match.arg(alternative)
    logic.operator <- match.arg(logic.operator)
    logic.operator <- ifelse(logic.operator=="or", "|", "&")
    predictions[which(is.nan(predictions))] <- NA
    observations[which(is.nan(observations))] <- NA
    cc.ix <- complete.cases(predictions, observations)
    predictions <- predictions[which(cc.ix)]
    observations <- observations[which(cc.ix)]
    N <- length(which(cc.ix))
    c.temp <- d <- u <- matrix(0, nrow = 1, ncol = N)
    c.d.seq <- NULL
    
    for (i in seq(from = 1, to = N - 1)) {
        for (j in seq(from = i + 1, to = N)) {
            pair <- c(i, j)
            iff <- abs(predictions[i] - predictions[j]) >= delta.pred | abs(observations[i] - observations[j]) >= delta.obs
            if(iff){ #add flag to replace 'or' behaviour with 'xor' behaviour
                pp <- (predictions[i] < predictions[j])
                oo <- (observations[i] < observations[j])
                if (pp == oo) {
                    c.temp[pair] <- c.temp[pair] + 1
                } else {
                    d[pair] <- d[pair] + 1
                }
                
            } else {
                if(outx){
                }else{
                    d[pair] <- d[pair] + 1
                }
            }
        }
    }
    
    C <- sum(c.temp)
    D <- sum(d)
    
    if (N < 3 || (C == 0 && D == 0)) {
        return(NA)
    }
    
    cindex <- C / (C + D)
    
    return(cindex)
}