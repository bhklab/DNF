rm(list=ls())

fake.sens.features <- matrix(rnorm(1000), nrow=100, ncol=10)
rownames(fake.sens.features) <- paste("feature", 1:nrow(fake.sens.features), sep="")
colnames(fake.sens.features) <- paste("drug", 1:ncol(fake.sens.features), sep="")
colnames(fake.sens.features) <- paste("originaldata", colnames(fake.sens.features), sep=":::")


fake.sens.similarities <- cor(fake.sens.features)
dimnames(fake.sens.similarities) <- list(paste("drug", 1:ncol(fake.sens.features), sep=""), paste("drug", 1:ncol(fake.sens.features), sep=""))

test_that("Adding Existing Drug", {
    new.compound.name <- "drug1"
    new.compound.features <- rnorm(100)
    
    expected.result <- list(updated.features=fake.sens.features, updated.similarities=fake.sens.similarities)
    
    expect_equal(AddNewSensitivityFeatures(fake.sens.features, fake.sens.similarities, new.compound.name,
                                           new.compound.features), expected.result)
})

test_that("Add New Drug With All Features in Common", {
    new.compound.name <- "drug11"
    new.compound.features <- rnorm(100)
    names(new.compound.features) <- rownames(fake.sens.features)
    
    expected.features <- cbind(fake.sens.features, new.compound.features)
    colnames(expected.features)[ncol(expected.features)] <- new.compound.name
    
    expected.similarities <- cor(expected.features, use = "pairwise.complete.obs")
    dimnames(expected.similarities) <- list(paste("drug", 1:ncol(expected.similarities), sep=""), paste("drug", 1:ncol(expected.similarities), sep=""))
    expected.similarities <- expected.similarities[order(rownames(expected.similarities)), 
                                                   order(colnames(expected.similarities))]
    
    colnames(expected.features)[colnames(expected.features) == new.compound.name] <- paste("new.dataset", new.compound.name, sep=":::")
    
    expected.result <- list(updated.features=expected.features, updated.similarities=expected.similarities)
    expect_equal(AddNewSensitivityFeatures(fake.sens.features, fake.sens.similarities, new.compound.name,
                                           new.compound.features), expected.result)
})

test_that("Add New Drug With Some Features in Common", {
    new.compound.name <- "drug11"
    new.compound.features <- rnorm(100)
    names(new.compound.features) <- rownames(fake.sens.features)
    uncommon.features <- sample(1:length(new.compound.features), 50)
    names(new.compound.features)[uncommon.features] <- paste("uncommonfeature", uncommon.features, sep="")
    
    expected.features <- cbind(fake.sens.features, new.compound.features)
    colnames(expected.features)[ncol(expected.features)] <- new.compound.name
    
    expected.features[uncommon.features, new.compound.name] <- NA
    
    expected.similarities <- cor(expected.features, use = "pairwise.complete.obs")
    dimnames(expected.similarities) <- list(paste("drug", 1:ncol(expected.similarities), sep=""), paste("drug", 1:ncol(expected.similarities), sep=""))
    expected.similarities <- expected.similarities[order(rownames(expected.similarities)), 
                                                   order(colnames(expected.similarities))]
    
    colnames(expected.features)[colnames(expected.features) == new.compound.name] <- paste("new.dataset", new.compound.name, sep=":::")
    
    expected.result <- list(updated.features=expected.features, updated.similarities=expected.similarities)
    expect_equal(AddNewSensitivityFeatures(fake.sens.features, fake.sens.similarities, new.compound.name,
                                           new.compound.features), expected.result)
})

test_that("Add New Drug with No Features in Common", {
    new.compound.name <- "drug11"
    new.compound.features <- rnorm(100)
    names(new.compound.features) <- paste("new.features", 1:100, sep="")
    
    expected.features <- fake.sens.features
    expected.similarities <- fake.sens.similarities
    expected.result <- list(updated.features=expected.features, updated.similarities=expected.similarities)
    
    expect_equal(AddNewSensitivityFeatures(fake.sens.features, fake.sens.similarities, new.compound.name, 
                                           new.compound.features), expected.result)
})