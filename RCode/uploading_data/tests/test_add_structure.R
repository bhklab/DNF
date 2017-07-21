rm(list=ls())

set.seed(123)
indices <- sample(1:nrow(lincs.meta), 10)

molecules <- rcdk::parse.smiles(lincs.meta$canonical_smiles[indices])
fingerprints <- lapply(molecules, rcdk::get.fingerprint, type="extended")
names(fingerprints) <- lincs.meta$pert_iname[indices]

fake.structure.features <- fingerprints
fake.structure.similarities <- fingerprint::fp.sim.matrix(fingerprints)
rownames(fake.structure.similarities) <- names(fake.structure.features)
colnames(fake.structure.similarities) <- names(fake.structure.features)

test_that("Adding Existing Drug", {
    new.compound.name <- names(fake.structure.features)[1]
    new.compound.features <- "C01[]C"
    
    expected.result <- list(updated.features=fake.structure.features, updated.similarities=fake.structure.similarities)
    
    expect_equal(AddNewStructureFeatures(fake.structure.features, fake.structure.similarities, new.compound.name,
                                            new.compound.features), expected.result)
})

test_that("Add New Drug With All Features in Common", {
    new.compound.name <- "drug11"
    new.compound.features <- lincs.meta$canonical_smiles[3000]
    
    new.compound.fingerprint <- rcdk::get.fingerprint(rcdk::parse.smiles(new.compound.features)[[1]], type="extended")
    
    expected.features <- fake.structure.features
    expected.features[[new.compound.name]] <- new.compound.fingerprint
    
    expected.similarities <- fingerprint::fp.sim.matrix(expected.features)
    rownames(expected.similarities) <- names(expected.features)
    colnames(expected.similarities) <- names(expected.features)
    
    expected.similarities <- expected.similarities[order(rownames(expected.similarities)), 
                                                   order(colnames(expected.similarities))]
    expected.result <- list(updated.features=expected.features, updated.similarities=expected.similarities)
    
    expect_equal(AddNewStructureFeatures(fake.structure.features, fake.structure.similarities, new.compound.name,
                                            new.compound.features), expected.result)
})