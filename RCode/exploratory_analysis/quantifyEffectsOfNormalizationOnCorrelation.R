rm(list=ls())

library(preprocessCore)

temp <- matrix(rexp(5000), nrow=100, ncol=50)
temp.quantile.normalized <- normalize.quantiles(temp)
temp.scaled <- scale(temp)

old.par <- par(mfrow=c(1,2))
colRamp <- colorRampPalette(c(3, "white", 2))(5)

plot(density(temp[,1]), col=colRamp[1], ylim=c(0,1))
for (i in 2:5) {lines(density(temp[, i]), col=colRamp[i])}

plot(density(temp.quantile.normalized[,1]), col=colRamp[1], ylim=c(0,1))
for (i in 2:5) {lines(density(temp.quantile.normalized[, i]), col=colRamp[i])}

cor(temp)
cor(temp.quantile.normalized)

pert.matrix <- matrix(0, nrow=100, ncol=50)
pert.indices <- sample(1:ncol(temp), 10)
pert.matrix[, pert.indices] <- pert.matrix[, pert.indices] + pert.indices
temp.perturbed <- temp + pert.matrix
temp.perturbed.quantile.normalized <- normalize.quantiles(temp.perturbed)

old.par <- par(mfrow=c(1,2))
colRamp <- colorRampPalette(c(3, "white", 2))(5)

plot(density(temp.perturbed[,1]), col=colRamp[1], ylim=c(0,1))
for (i in 2:5) {lines(density(temp.perturbed[, i]), col=colRamp[i])}

plot(density(temp.perturbed.quantile.normalized[,1]), col=colRamp[1], ylim=c(0,1))
for (i in 2:5) {lines(density(temp.perturbed.quantile.normalized[, i]), col=colRamp[i])}

cor(temp.perturbed)
cor(temp.perturbed.quantile.normalized)

aff.perturbed <- SNFtool::affinityMatrix(1 - cor(temp.perturbed))
aff.perturbed.quantile.normalized <- SNFtool::affinityMatrix(1 - cor(temp.perturbed.quantile.normalized))
correlations <- vapply(1:ncol(aff.perturbed), function(i, unnormalized, normalized) {
    correlation <- cor(unnormalized[, i], normalized[, i])
    
    correlation
}, unnormalized=aff.perturbed, normalized=aff.perturbed.quantile.normalized, FUN.VALUE = numeric(1))

load("PSets/CTRPv2.RData")
aucs <- summarizeSensitivityProfiles(pSet = CTRPv2, sensitivity.measure = "auc_recomputed")
aucs.quantile.normalized <- normalize.quantiles(aucs)

cor.aucs <- cor(aucs, use = "pairwise.complete.obs")
cor.aucs.quantile.normalized <- cor(aucs.quantile.normalized, use = "pairwise.complete.obs")
cor.aucs[is.na(cor.aucs)] <- 0
cor.aucs.quantile.normalized[is.na(cor.aucs.quantile.normalized)] <- 0

aff.aucs <- SNFtool::affinityMatrix(1 - cor.aucs)
aff.aucs.quantile.normalized <- SNFtool::affinityMatrix(1 - cor.aucs.quantile.normalized)

correlations <- vapply(1:ncol(aff.aucs), function(i, unnormalized, normalized) {
    correlation <- cor(unnormalized[, i], normalized[, i], method = "spearman")
    
    correlation
}, unnormalized=aff.aucs, normalized=aff.aucs.quantile.normalized, FUN.VALUE = numeric(1))
names(correlations) <- colnames(aff.aucs)

