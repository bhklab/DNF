imaging.data <- read.delim("Data/Broad.HG005032.ProfilingData/imaging/cdrp.imaging.profiles.txt", stringsAsFactors = FALSE)

random.sample <- sample(1:ncol(imaging.data), 30, replace=FALSE)
boxplot(imaging.data[, random.sample])

colRamp <- colorRampPalette(c(3, "white", 2))(20)
plot(density(imaging.data[, 2]), lwd=3, col=colRamp[1])

imaging.data.scaled <- scale(imaging.data[, random.sample])
imaging.data.scaled <- normalize.quantiles(as.matrix(imaging.data[, random.sample]))

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
