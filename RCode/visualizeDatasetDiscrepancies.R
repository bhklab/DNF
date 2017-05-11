visualize_dataset_discrepancies <- function(dataset_pairs, unlisted_num_samples) {
    d <- data.frame(x = unlist(dataset_pairs), 
                    grp = rep(names(dataset_pairs), times = sapply(dataset_pairs, length)),
                    n_sample = unlisted_num_samples)
    p <- ggplot(d,aes(x = grp, y = x)) + geom_boxplot()
    # p + geom_jitter(aes(colour = n_sample)) + scale_colour_gradientn(colours = terrain.colors(9)) + ggtitle('Correlation of same Drug With Itself Across Dataset Pairs')
    br <- seq(0, 900, length.out = 10)
    res <- p + geom_jitter(aes(colour = n_sample)) + scale_colour_gradientn(colours = rainbow(3), 
                         na.value = "white", 
                         breaks = br) + ggtitle('Correlation of same Drug With Itself Across Dataset Pairs') + labs(x = "Dataset Pair", y = "Correlation")
    print(res)
}