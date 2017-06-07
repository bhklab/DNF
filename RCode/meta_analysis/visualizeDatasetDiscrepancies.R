VisualizeDatasetDiscrepancies <- function(dataset.pairs, unlisted.num.samples) {
    # Creates a series of boxplots, where each boxplot corresponds to correlations
    # between a pair of datasets. For example, dataset 1 and dataset 2, the boxplot 
    # would contain correlations between each drug x in dataset 1, and the same 
    # drug x in dataset 2.
    #
    # Args: 
    #   dataset.pairs: A list of dataset pairs where the value of each list element
    #                  is a vector of correlations. Each correlation is between a drug
    #                  and itself between those dataset pairs.
    #   unlisted.num.samples: A vector containing the number of samples used to compute
    #                         the various correlations that are found in dataset.pairs.
    #
    d <- data.frame(x = unlist(dataset.pairs), 
                    grp = rep(names(dataset.pairs), times = sapply(dataset.pairs, length)),
                    n.sample = unlisted.num.samples)
    p <- ggplot(d,aes(x = grp, y = x)) + geom_boxplot() + scale_x_discrete(labels=levels(d$grp))
    # p + geom_jitter(aes(colour = n_sample)) + scale_colour_gradientn(colours = terrain.colors(9)) + ggtitle('Correlation of same Drug With Itself Across Dataset Pairs')
    br <- seq(0, 900, length.out = 10)
    res <- p + geom_jitter(aes(colour = n.sample)) + scale_colour_gradientn(colours = rainbow(3), 
                         na.value = "white", 
                         breaks = br) + ggtitle('Correlation of same Drug With Itself Across Dataset Pairs') + labs(x = "Dataset Pair", y = "Correlation")
    print(res)
}