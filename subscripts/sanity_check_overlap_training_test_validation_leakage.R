setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v2.rds")

# Sanity check overlap between sets (CV leakage) ----
meta[set %in% c("test", "validation", "training"), {
  # Import regions
  .c <- .SD[, {
    importBed(bed_file)
  }, .(set, bed_file)]
  # Split sets
  .c <- split(.c, .c$set)
  # Overlaps between sets
  for(i in seq(length(.c)-1)) {
    for(j in (i+1):length(.c)) {
      if(nrow(intersectBed(.c[[i]], .c[[j]])))
        stop(paste(paste0(unlist(.BY), collapse= "+"), "-> Overlap were found between sets"))
    }
  }
  # Progress
  print(paste0(.GRP, "/", .NGRP))
}, .(ID, dataset, tissue, fold)]
