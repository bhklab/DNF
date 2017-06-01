badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"

# Clean up badchars in names and capitalize
manual.mapping <- read.csv("Data/mapping_manual_curation_inchikeys.csv")
manual.mapping$drugnames <- toupper(manual.mapping$drugnames)
manual.mapping$drugnames <- gsub(badchars, "", manual.mapping$drugnames)
manual.mapping$pert_iname <- toupper(manual.mapping$pert_iname)
manual.mapping$pert_iname <- gsub(badchars, "", manual.mapping$pert_iname)

# Clean up badchars in names and capitalize
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)
lincs.meta$pert_iname <- toupper(lincs.meta$pert_iname)
lincs.meta$pert_iname <- gsub(badchars, "", lincs.meta$pert_iname)

sens.data <- readRDS("Data/combined_sens_adjusted_diag_datasets_with.RData")

# See intersection between sensitivity data and manual mapping drug names
intersection.between.sens.manual.drugs <- intersect(rownames(sens.data), manual.mapping$drugnames)
intersection.between.sens.manual.iname <- intersect(rownames(sens.data), manual.mapping$pert_iname)

# See intersection between lincs pert_iname and sensitivity data
intersection.between.sens.lincs <- intersect(rownames(sens.data), lincs.meta$pert_iname)

# Combine all the drugnames and pert_iname of the drugs that intersected between our sensitivity data and the file
# that Nehme sent us.
combined.intersection <- union(intersection.between.sens.manual.drugs, intersection.between.sens.manual.iname)

# Take the difference between the combined.intersection drugs and the intersection between sensitivity and lincs.
# This gives us the drugs that we will gain from using the new matching method. The intersection between 
# sens and lincs is 321 drugs. The combined.intersection drugs (the ones between sensitivity data and 
# the file that Nehme sent) is 302 drugs. Nehme's file contains only drugs that are found in L1000 as well,
# so taking the difference between combined.intersection and intersection.between.sens.lincs will give us
# the additional drugs we will get by incorporating Nehme's method.
setdiff(combined.intersection, intersection.between.sens.lincs)

complete.intersection.between.sens.manual <- unique(union(intersection.between.sens.manual.drugs,
                                                          intersection.between.sens.manual.iname))

# See intersection between lincs pert_iname and manual mapping pert_iname
intersection.between.lincs.manual.iname <- intersect(lincs.meta$pert_iname, manual.mapping$pert_iname)

# Parse out only first part of inchi key from lincs meta. Inchi keys have form xx_yy_zz,
# and we care about xx in this case
splitted <- strsplit(lincs.meta$inchi_key, split="=")
new.inchi <- sapply(splitted, function(x) {x[2]})
splitted <- strsplit(new.inchi, split="-")
new.inchi <- sapply(splitted, function(x) {x[1]})
lincs.meta$inchi_key <- new.inchi

# See intersection between inchi keys from manual mapping and those in LINCS meta
intersection.between.ichi.manual <- intersect(new.inchi, manual.mapping$inchi_key)


lincs.meta.subsetted <- lincs.meta[lincs.meta$inchi_key %in% manual.mapping$inchi_key,]

setdiff(lincs.meta.subsetted$pert_iname, rownames(sens.data))
