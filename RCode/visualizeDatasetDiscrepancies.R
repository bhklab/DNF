VisualizeDatasetDiscrepancies <- function(dataset.pairs, unlisted.num.samples) {
    d <- data.frame(x = unlist(dataset.pairs), 
                    grp = rep(names(dataset.pairs), times = sapply(dataset.pairs, length)),
                    n.sample = unlisted.num.samples)
    p <- ggplot(d,aes(x = grp, y = x)) + geom_boxplot()
    # p + geom_jitter(aes(colour = n_sample)) + scale_colour_gradientn(colours = terrain.colors(9)) + ggtitle('Correlation of same Drug With Itself Across Dataset Pairs')
    br <- seq(0, 900, length.out = 10)
    res <- p + geom_jitter(aes(colour = n.sample)) + scale_colour_gradientn(colours = rainbow(3), 
                         na.value = "white", 
                         breaks = br) + ggtitle('Correlation of same Drug With Itself Across Dataset Pairs') + labs(x = "Dataset Pair", y = "Correlation")
    print(res)
}