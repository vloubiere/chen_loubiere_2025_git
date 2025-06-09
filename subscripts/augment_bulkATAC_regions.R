setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
source("git_deepATAC/function/augmentation_function_tiling_sliding_window.R")
require(vlfunctions)

# Import folds ----
dat <- readRDS("db/folds/bulkATAC_folds.rds")

# Import bw metadata for coverage ----
bw <- as.data.table(readxl::read_xlsx("Rdata/metadata_ATACSeq.xlsx"))[dataset=="bulkENCODE"]

# VISTA TILES -------------------------------------------------------------------#
vista <- dat[group=="vista"]
vista[, c("seqnames", "start", "end"):= tstrsplit(coor, ":|-", type.convert = T)] # Retrieve mouse/human coordinates.
augVISTA <- augSlideFunction(vista, width = 1001, step= 50)
# Extract sequences
augVISTA[genome=="hg38", seq:= vl_getSequence(.SD, "hg38")]
augVISTA[genome=="mm10", seq:= vl_getSequence(.SD, "mm10")]
augVISTA[, ID:= paste0(peakID, ":", strand, "_", shift)]

# GLOBALLY OPEN -----------------------------------------------------------------#
GO <- dat[group=="ATAC" & class=="globallyOpen"]
# Augment 10x
augGO <- augTile(GO, width = 1001, shifts = seq(-400, 400, 200))
# Extract coverage and sequences
augGO[, paste0(bw$tissue, ".cov"):= lapply(bw$bw, function(x) log2(vl_bw_coverage(.SD, x)+1))]
augGO[, seq:= vl_getSequence(augGO, "mm10")]
augGO[, ID:= paste0(peakID, ":", strand, "_", shift)]

# TISSUE-SPECIFIC ---------------------------------------------------------------#
TS <- dat[group=="ATAC" & class!="globallyOpen"]
# Augment 34x
augTS <- augTile(TS, width = 1001, shifts = seq(-400, 400, 50))
# Extract coverage and sequences
augTS[, paste0(bw$tissue, ".cov"):= lapply(bw$bw, function(x) log2(vl_bw_coverage(.SD, x)+1))]
augTS[, seq:= vl_getSequence(augTS, "mm10")]
augTS[, ID:= paste0(peakID, ":", strand, "_", shift)]

# NEGATIVE CONTROLS --------------------------------------------------------------#
CTL <- dat[group=="ctl"]
# No augmentation
CTL[, shift:= 0L]
set.seed(1)
CTL[, strand:= sample(c("+", "-"), .N, replace = T)]
CTL <- vl_resizeBed(CTL, "center", 500, 500)
# Extract coverage and sequences (no augmentation)
CTL[, paste0(bw$tissue, ".cov"):= lapply(bw$bw, function(x) log2(vl_bw_coverage(.SD, x)+1))]
CTL[, seq:= vl_getSequence(CTL, "mm10")]
CTL[, ID:= paste0(peakID, ":", strand, "_", shift)] # shift equals 0

# Subsample before augmenting 10x ------------------------------------------------#
set.seed(1)
augCTL <- CTL[sample(.N, 1e5)]
augCTL$strand <- NULL
augCTL <- augTile(augCTL, width = 1001, shifts = seq(-400, 400, 200))
# Extract coverage and sequences
augCTL[, paste0(bw$tissue, ".cov"):= lapply(bw$bw, function(x) log2(vl_bw_coverage(.SD, x)+1))]
augCTL[, seq:= vl_getSequence(augCTL, "mm10")]
augCTL[, ID:= paste0(peakID, ":", strand, "_", shift)]

# Remove repeated sequences ----
# Here, the order matters, as only the first (top) instance of duplicated sequences will be kept
cmb <- list("db/augmentation/bulkATAC/vistaTiles_1001_50.rds"= augVISTA,
            "db/augmentation/bulkATAC/tissueSpecific_34x.rds"= augTS,
            "db/augmentation/bulkATAC/globallyOpen_10x.rds"= augGO,
            "db/augmentation/bulkATAC/globallyClosed_noAUG.rds"= CTL,
            "db/augmentation/bulkATAC/globallyClosed_10x.rds"= augCTL)
cmb <- rbindlist(cmb, idcol = "file", fill = T)
# Remove long NNNN stretches ----
cmb[, longN:= grepl(paste0(rep("N", 250), collapse= ""), cmb$seq)]
# Remove duplicated sequences ----
cmb[, mult:= uniqueN(ID)>1, seq]
cmb[(mult), mult:= ID!=ID[1], seq] # Can be kept once 
# Report and remove ----
cmb[, .(.N,
        nIDs= length(unique(peakID)),
        longN= sum(longN),
        longNpeakIDs= length(unique(peakID[longN])),
        repeated= sum(mult),
        repeatedIDs= length(unique(peakID[mult]))), file]
cmb <- cmb[(!longN & !mult), !c("longN", "mult")]
# Save ----
cmb[, {
  # Simplify (remove empty columns)
  .c <- data.table::copy(.SD)
  empty <- sapply(.c, function(x) all(is.na(x)))
  # Save
  saveRDS(.c[, !empty, with= F],
          file)
  print(file)
}, file]