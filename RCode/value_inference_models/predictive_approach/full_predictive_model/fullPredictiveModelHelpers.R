CreateModels <- function(input.layers, output.layer, prediction.matrix, model.name) {
    models <- list()
    
    #feature.names <- colnames(input.layers[['sensitivity']])
    if (length(input.layers) > 1) {
        feature.names <- Reduce(function(x, y) 
            {union(colnames(x), colnames(y))}, input.layers)
        feature.names <- sapply(feature.names, function(x) {paste('`', x, '`', sep="")})
    } else {
        feature.names <- colnames(input.layers[[1]])
        feature.names <- sapply(feature.names, function(x) {paste('`', x, '`', sep="")})
    }
    
    
    # Build a regression model for each image feature and store it in the models list
    for (i in 1:ncol(output.layer)) {
        target.name <- colnames(output.layer)[i]
        
        target.name.temp <- paste(target.name, " ~ ", sep="")
        
        prediction.formula <- as.formula(paste(target.name.temp, paste(feature.names, sep="", collapse = " + "), sep=""))
        
        ctrl <- trainControl(method="none")
        models[[target.name]] <- train(prediction.formula, data=prediction.matrix, method=model.name,
                                       na.action=na.pass, trControl = ctrl)
        print(i)
    }    
    
    models
}

MakePredictions <- function(models, drugs.to.predict, feature.names, prediction.matrix) {
    predicted.results <- matrix(0, nrow=length(drugs.to.predict), ncol=length(feature.names))
    rownames(predicted.results) <- drugs.to.predict
    colnames(predicted.results) <- feature.names
    
    for (i in 1:length(models)) {
        feature.name <- names(models)[i]
        predicted.results[, feature.name] <- predict(models[[feature.name]], prediction.matrix)
    }
    
    predicted.results
}

CreateSensitivityFeatures <- function(datasets, badchars, names.to.intersect) {
    datasets <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)
    
    cell.lines <- sort(unique(unlist(sapply(datasets, function(x) { return (cellNames(x)) }))))
    drugs <- unique(unlist(sapply(datasets, function(x) { return (paste(pSetName(x), gsub(badchars, "",toupper(drugNames(x))), sep=":::")) })))
    
    ### Merge datasets into one large matrix where rows are drugs (prepended with the corresponding dataset),
    ### and columns are cell lines
    aucs.all <- CreateAucsAll(datasets, drugs, cell.lines = cell.lines, badchars = badchars)
    
    splitted <- strsplit(rownames(aucs.all), split=":::")
    sens.names <- sapply(splitted, function(x) x[2])
    overlap <- sort(Reduce(intersect, list(sens.names, names.to.intersect)))
    
    # Create placeholder for final AUC matrix
    aucs.all.final <- matrix(0, nrow=length(overlap), ncol=ncol(aucs.all))
    rownames(aucs.all.final) <- overlap
    colnames(aucs.all.final) <- colnames(aucs.all)
    
    # Combine AUCs for drugs found in multiple datasets by averaging values. Also, only take subset of columns
    # having less than 20% NA values
    aucs.all.final <- AverageAUCS(overlap, aucs.all, aucs.all.final)
    columns.to.keep <- colSums(is.na(aucs.all.final)) < (0.2 * nrow(aucs.all.final))
    aucs.all.final <- aucs.all.final[, columns.to.keep]
    
    # Impute remaining missing value
    c1 <- makeCluster(8)
    registerDoParallel(c1)
    aucs.all.final <- missForest(aucs.all.final, maxiter = 10, parallelize = "variables")$ximp
    stopCluster(c1)
    
    list(sens.data.final=aucs.all.final, overlap=overlap)
}

CreateStructureFeatures <- function(lincs.meta) {
    fingerprints <- StructureDataFlexible(lincs.meta)
    
    c1 <- makeCluster(8)
    registerDoParallel(c1)
    
    res <- foreach(i=1:length(fingerprints)) %dopar% {
        placeholder <- numeric(1024)
        fp <- fingerprints[[i]]
        bits <- fp@bits
        
        placeholder[bits] <- 1
        placeholder
    }
    
    stopCluster(c1)
    
    df.fingerprints <- matrix(unlist(res), nrow=length(res), byrow=T)
    rownames(df.fingerprints) <- names(fingerprints)
    
    df.fingerprints
}

CreatePerturbationFeatures <- function(pert.file.name, badchars) {
    df.pert <- PerturbationDataFlexible(pert.file.name, lincs.meta)
    df.pert <- t(df.pert)
    rownames(df.pert) <- toupper(rownames(df.pert))
    rownames(df.pert) <- gsub(badchars, "", rownames(df.pert))
    
    df.pert <- df.pert[rowSums(is.na(df.pert)) < (0.9 * ncol(df.pert)), ]
    colnames(df.pert) <- gsub("-", "_", colnames(df.pert))
    
    df.pert
}

