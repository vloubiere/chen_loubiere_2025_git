setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
dat <- readRDS("Rdata/metadata_ATACSeq_models_processed.rds")

# Import mm10 cordinates VISTA ----
coor <- readRDS("db/peaks/vista_tiles_clean.rds")

# Check that there is no overlap between test and training/valid for each fold
dat[tissue==tissue[1] & ID==ID[1], # Folds are shared between tissues and models, no need for individual testing
{
  # Final test set transfer learning
  testVISTA <- as.data.table(rtracklayer::import(bed_VISTA_test))
  testVISTA[grepl("^hs", name), ID:= tstrsplit(name, ":", keep= 1)]
  testVISTA[coor, c("seqnames", "start", "end"):= .(i.seqnames, i.start, i.end), on= "ID==peakID"] # Retrieve mm10 coor for vista hs sequences
  
  # Valid/test sets
  validVISTA <- as.data.table(rtracklayer::import(bed_VISTA_validation))
  trainVISTA <- as.data.table(rtracklayer::import(bed_VISTA_training))
  validATAC <- as.data.table(rtracklayer::import(bed_ATAC_validation))
  trainATAC <- as.data.table(rtracklayer::import(bed_ATAC_training))
  others <- rbind(validVISTA,
                  trainVISTA,
                  validATAC,
                  trainATAC)
  others[grepl("^hs", name), ID:= tstrsplit(name, ":", keep= 1)]
  others[coor, c("seqnames", "start", "end"):= .(i.seqnames, i.start, i.end), on= "ID==peakID"] # Retrieve mm10 coor for vista hs sequences
  
  # Overlaps between VISTA test and rest (except ATAC test)
  VISTA <- vl_intersectBed(testVISTA[!is.na(seqnames)],
                           others[!is.na(seqnames)])
  print(paste(nrow(VISTA), "overlaps found"))
}, .(fold,
     bed_ATAC_test, bed_ATAC_validation, bed_ATAC_training,
     bed_VISTA_test, bed_VISTA_validation, bed_VISTA_training)]
