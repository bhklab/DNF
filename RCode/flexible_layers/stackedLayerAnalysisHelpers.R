CreateBaseFileName <- function() {
    base.dir <- "Output/auc_p_flex"
    
    if(!file.exists(base.dir)) {
        dir.create(base.dir)
    }
    
    file.name <- paste(base.dir, "/", sep="")
}

AddLayersToFileName <- function(file.name, sensitivity.file.name="", use.sensitivity=FALSE,
                                use.perturbation=FALSE, use.structure=FALSE, use.luminex=FALSE,
                                use.imaging=FALSE) {
    if (use.sensitivity) {
        file.name <- paste(file.name, "sens", sep="")
        sensitivity.file.name <- strsplit(sensitivity.file.name, "/")[[1]][2]
        sensitivity.file.name <- strsplit(sensitivity.file.name, "[.]")[[1]][1]        
        file.name <- paste(file.name, sensitivity.file.name, sep="_")
    }
    
    if (use.perturbation) {
        file.name <- paste(file.name, "pert", sep="_")
    }
    
    if (use.structure) {
        file.name <- paste(file.name, "strc", sep="_")
    }
    
    if (use.luminex) {
        file.name <- paste(file.name, "luminex", sep="_")
    }
    
    if (use.imaging) {
        file.name <- paste(file.name, "imaging", sep="_")
    }
    
    file.name
}

AddDrugTargetBenchmarksToFileName <- function(file.name, use.ctrpv2, use.clue, use.chembl,
                                              use.dbank, use.dtc) {
    if (use.ctrpv2) {
        file.name <- paste(file.name, "ctrpv2", sep="_")
    }
    
    if (use.clue) {
        file.name <- paste(file.name, "clue", sep="_")
    }
    
    if (use.chembl) {
        file.name <- paste(file.name, "chembl", sep="_")
    }
    
    if (use.dbank) {
        file.name <- paste(file.name, "dbank", sep="_")
    }
    
    if (use.dtc) {
        file.name <- paste(file.name, "dtc", sep="_")
    }
    
    file.name
}

CreateTargetROCFileName <- function(sensitivity.file.name="", use.sensitivity=FALSE,
                              use.perturbation=FALSE, use.structure=FALSE, use.luminex=FALSE,
                              use.imaging=FALSE, use.ctrpv2=FALSE, use.clue=FALSE,
                              use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE) {
    file.name <- CreateBaseFileName()
    file.name <- AddLayersToFileName(file.name=file.name, sensitivity.file.name=sensitivity.file.name,
                                     use.sensitivity=use.sensitivity, use.perturbation=use.perturbation,
                                     use.structure=use.structure, use.luminex=use.luminex,
                                     use.imaging=use.imaging)
    
    file.name <- paste(file.name, "target", sep="_")
    file.name <- AddDrugTargetBenchmarksToFileName(file.name, use.ctrpv2=use.ctrpv2,
                                                   use.clue=use.clue, use.chembl=use.chembl,
                                                   use.dbank=use.dbank, use.dtc=use.dtc)
    
    file.name <- paste(file.name, "pdf", sep=".")
}

CreateATCROCFileName <- function(sensitivity.file.name="", atc.benchmark.name="",
                                 use.sensitivity=FALSE,
                                 use.perturbation=FALSE, use.structure=FALSE, use.luminex=FALSE,
                                 use.imaging=FALSE) {
    file.name <- CreateBaseFileName()
    
    file.name <- AddLayersToFileName(file.name=file.name, sensitivity.file.name=sensitivity.file.name,
                                     use.sensitivity=use.sensitivity, use.perturbation=use.perturbation,
                                     use.structure=use.structure, use.luminex=use.luminex,
                                     use.imaging=use.imaging)
    
    file.name <- paste(file.name, "atc", sep="_")
    file.name <- paste(file.name, atc.benchmark.name, sep="_")
    
    file.name <- paste(file.name, "pdf", sep=".")
}

CreateGMTFileName <- function(use.sensitivity=FALSE,
                              use.perturbation=FALSE, use.structure=FALSE, use.luminex=FALSE,
                              use.imaging=FALSE, use.ctrpv2=FALSE, use.clue=FALSE,
                              use.chembl=FALSE, use.dbank=FALSE, use.dtc=FALSE) {
    file.name <- CreateBaseFileName()
    
    file.name <- paste(file.name, "communities", sep="")
    
    file.name <- AddLayersToFileName(file.name=file.name, sensitivity.file.name="",
                                     use.sensitivity=use.sensitivity, use.perturbation=use.perturbation,
                                     use.structure=use.structure, use.luminex=use.luminex,
                                     use.imaging=use.imaging)
    
    file.name <- AddDrugTargetBenchmarksToFileName(file.name, use.ctrpv2=use.ctrpv2,
                                                   use.clue=use.clue, use.chembl=use.chembl,
                                                   use.dbank=use.dbank, use.dtc=use.dtc)
    
    file.name <- paste(file.name, "RData", sep=".")
}