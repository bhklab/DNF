lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
imaging.meta <- read.delim("Data/Broad.HG005032.ProfilingData/luminex/cdrp.l1000.meta.cpd.txt", 
                           stringsAsFactors = FALSE)

lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)
lincs.meta$pert_id <- toupper(lincs.meta$pert_id)
lincs.meta$pert_id <- gsub(badchars, "", lincs.meta$pert_id)


imaging.meta$name <- toupper(imaging.meta$name)
imaging.meta$name <- gsub(badchars, "", imaging.meta$name)

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"
# Get common drugs between imaging data and LINCs
pert.id.inter <- which(imaging.meta$name %in% lincs.meta$pert_id)
pert.name.inter <- which(imaging.meta$name %in% lincs.meta$pert_iname)
pert.name.inter <- imaging.meta$name[pert.name.inter]
pert.name.inter <- pert.name.inter[!duplicated(pert.name.inter)]
pert.name.inter <- match(pert.name.inter, imaging.meta$name)

combined <- unique(union(pert.id.inter, pert.name.inter))
temp.combined <- imaging.meta$name[combined]
temp.combined <- temp.combined[!duplicated(temp.combined)]
combined <- match(temp.combined, imaging.meta$name)

any(duplicated(imaging.meta$name))
sum(duplicated(imaging.meta$name))
any(duplicated(imaging.meta[combined, "name"]))

sens_data <- readRDS("Data/combined_sens_adjusted_diag_datasets_with.RData")
lincs.names <- lincs.meta[pert.name.inter, "pert_iname"]

lincs.names <- toupper(lincs.names)
lincs.names <- gsub(badchars, "", lincs.names)

sensitivity.inter <- intersect(rownames(sens_data), lincs.names)
