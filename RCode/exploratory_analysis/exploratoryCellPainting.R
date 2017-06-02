imaging.subsetted <- readRDS("Data/imaging_subsetted.RData")
imaging.subsetted <- t(imaging.subsetted)

random.sample <- sample(1:ncol(imaging.subsetted), 20, replace=FALSE)
boxplot(imaging.subsetted[, random.sample])

imaging.subsetted <- scale(imaging.subsetted)
imaging.subsetted <- imaging.subsetted[, colSums(is.na(imaging.subsetted)) != nrow(imaging.subsetted)]
princ <- prcomp(imaging.subsetted)
n.comp <- 10

df.components <- predict(princ, newdata=imaging.subsetted)[,1:n.comp]
boxplot(df.components)

plot(density(df.components[,1]), col=1)

for (i in 1:ncol(df.components)) {
    lines(density(df.components[,i]), col=i)
}

drugs.of.interest <- df.components[!startsWith(rownames(df.components), 'BRD'), ]
drugs.of.interest <- drugs.of.interest[sample(1:nrow(drugs.of.interest), 100, replace = FALSE),]

scatterplot3d(x=drugs.of.interest[,1], y=drugs.of.interest[,2], drugs.of.interest[,3])
stats::heatmap(drugs.of.interest, Rowv=NA, Colv=NA)
