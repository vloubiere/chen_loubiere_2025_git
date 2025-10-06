setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v2.rds")
meta <- meta[ID=="model1_bulkATAC_tsx3Aug_2xBal_noW"]
meta <- meta[(dataset=="accessibility" & set=="testCenterActBins") | (dataset=="activity" & set=="test")]

# Import bed files ----
bed <- meta[, importBed(bed_file), .(bed_file)]
bed[, genome:= ifelse(grepl("^hs", name), "hg38", "mm10")]
bed[, c("name", "label"):= tstrsplit(name, "__")]

# Collapse
bed <- bed[, collapseBed(.SD), genome]
bed[, strand:= "+"]
bed[, name:= paste0(seqnames, ":", start, "-", end, ":", strand, "__", genome)]
exportBed(bed[, !"genome"], "db/bed/collapsed_test_set_paper.bed")
