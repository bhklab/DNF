rm(list=ls())

fake.pert.features <- matrix(rnorm(1000), nrow=100, ncol=10)
rownames(fake.pert.features) <- paste("feature", 1:nrow(fake.pert.features), sep="")
colnames(fake.pert.features) <- paste("drug", 1:ncol(fake.pert.features), sep="")

fake.pert.similarities <- cor(fake.pert.features)

test_that("Adding Existing Drug", {
    new.compound.name <- "drug1"
    new.compound.features <- rnorm(100)
    
    expected.result <- list(updated.features=fake.pert.features, updated.similarities=fake.pert.similarities)
    
    expect_equal(AddNewPerturbationFeatures(fake.pert.features, fake.pert.similarities, new.compound.name,
                                           new.compound.features), expected.result)
})

test_that("Add New Drug With All Features in Common", {
    new.compound.name <- "drug11"
    new.compound.features <- rnorm(100)
    names(new.compound.features) <- rownames(fake.pert.features)
    
    expected.features <- cbind(fake.pert.features, new.compound.features)
    colnames(expected.features)[ncol(expected.features)] <- new.compound.name
    
    expected.similarities <- cor(expected.features, use = "pairwise.complete.obs")
    expected.similarities <- expected.similarities[order(rownames(expected.similarities)), 
                                                   order(colnames(expected.similarities))]
    expected.result <- list(updated.features=expected.features, updated.similarities=expected.similarities)

    expect_equal(AddNewPerturbationFeatures(fake.pert.features, fake.pert.similarities, new.compound.name,
                                           new.compound.features), expected.result)
})

test_that("Add New Drug With Some Features in Common", {
    new.compound.name <- "drug11"
    new.compound.features <- rnorm(100)
    names(new.compound.features) <- rownames(fake.pert.features)
    uncommon.features <- sample(1:length(new.compound.features), 50)
    names(new.compound.features)[uncommon.features] <- paste("uncommonfeature", uncommon.features, sep="")
    
    expected.features <- cbind(fake.pert.features, new.compound.features)
    colnames(expected.features)[ncol(expected.features)] <- new.compound.name
    
    expected.features[uncommon.features, new.compound.name] <- NA
    
    expected.similarities <- cor(expected.features, use = "pairwise.complete.obs")
    expected.similarities <- expected.similarities[order(rownames(expected.similarities)), 
                                                   order(colnames(expected.similarities))]
    
    expected.result <- list(updated.features=expected.features, updated.similarities=expected.similarities)

    expect_equal(AddNewPerturbationFeatures(fake.pert.features, fake.pert.similarities, new.compound.name,
                                           new.compound.features), expected.result)
})

test_that("Add New Drug with No Features in Common", {
    new.compound.name <- "drug11"
    new.compound.features <- rnorm(100)
    names(new.compound.features) <- paste("new.features", 1:100, sep="")
    
    expected.features <- fake.pert.features
    expected.similarities <- fake.pert.similarities
    expected.result <- list(updated.features=expected.features, updated.similarities=expected.similarities)
    
    expect_equal(AddNewPerturbationFeatures(fake.pert.features, fake.pert.similarities, new.compound.name, 
                                           new.compound.features), expected.result)
})
