source("RCode/uploading_data/uploadHelpers.R")

feature.file.paths <- list(sensitivity="Data/uploading_features/sensitivity/aucs_all.RData", 
                           perturbation="Data/uploading_features/perturbation/pert_features.RData",
                           structure="Data/uploading_features/structure/structure_features.RData",
                           imaging="")

similarity.file.paths <- list(sensitivity="Data/uploading_features/sensitivity/aucs_cor2.RData", 
                              perturbation="Data/uploading_features/perturbation/perturbation_similarities.RData",
                                structure="Data/uploading_features/structure/structure_similarities.RData",
                                imaging="")

new.compound.name <- "fake_drug"
fake.pert.features <- rnorm(978)
names(fake.pert.features) <- rownames(temp)

new.compound.features <- list(perturbation=fake.pert.features)

UploadFeaturesForNewCompound(new.compound.features, new.compound.name)

new.compound.name <- "fake_drug"
fake.sens.features <- rnorm(length(rownames(temp)))
names(fake.sens.features) <- rownames(temp)

new.compound.features <- list(sensitivity=fake.sens.features)

UploadFeaturesForNewCompound(new.compound.features, new.compound.name)

new.compound.name <- "fake_drug"
fake.strc.smiles <- "COCC(=O)N1CC[C@@]2(CC1)CN(Cc1ccccc1Cl)[C@@H](CO)c1[nH]c3cc(OC)ccc3c21"

new.compound.features <- list(structure=fake.strc.smiles)

UploadFeaturesForNewCompound(new.compound.features, new.compound.name)
