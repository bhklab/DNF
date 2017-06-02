load("PSets/CCLE_hs.RData")
load("PSets/gCSI_hs.RData")
load("PSets/CTRPv2.RData")
load("PSets/GDSC1000.RData")
load("PSets/FIMM.RData")
badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[_]|[ ]"


datasets <- list(CCLE, gCSI, GDSC1000, CTRPv2, FIMM)

drug.names <- unique(unlist(sapply(datasets, function(x) { return (gsub(badchars, "",toupper(drugNames(x)))) })))
drug.file <- read.csv("Data/drugs_with_ids.csv", stringsAsFactors = FALSE)

drug.file$unique.drugid <- toupper(drug.file$unique.drugid)
drug.file$unique.drugid <- gsub(badchars, "", drug.file$unique.drugid)

intersection <- drug.file[drug.file$unique.drugid %in% drug.names, ]
sum(is.na(intersection$inchikey))
lincs.meta <- read.csv("Data/LINCS.csv", stringsAsFactors = FALSE)

any(is.na(lincs.meta$inchi_key))
View(drugInfo(CCLE))
