setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
source("git_deepATAC/function/augmentation_function_tiling_sliding_window.R")
require(vlfunctions)

# Import ATAC-Seq peaks, vista tiles, control regions and compute overlaps ----
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
vista[, start:= start-100]# Extended to reflect later augmentation
vista[, end:= end+100]# Extended to reflect later augmentation
ATAC <- readRDS("db/peaks/bulkENCODE_confident_ATAC_peaks.rds")
ATAC <- vl_resizeBed(ATAC, "center", 1400, 1400, genome = "mm10") # Extended to reflect later augmentation
ctl <- readRDS("db/peaks/bulkENCODE_control_regions.rds")
ctl <- vl_resizeBed(ctl, "center", 1400, 1400, genome = "mm10") # Extended to reflect later augmentation

# Combine and compute overlaps ----
cmb <- list(vista= vista,
            ATAC= ATAC,
            ctl= ctl)
cmb <- rbindlist(cmb, fill = T, idcol= "group")
cmb[!is.na(seqnames), ov:= vl_collapseBed(cmb[!is.na(seqnames)], return.idx.only = T)]
cmb[is.na(ov), ov:= max(cmb$ov, na.rm = T)+seq(.N)] # Human-specific sequences have no mm10 coordinates
if(nrow(cmb) != length(unique(cmb$peakID)))
  stop("peakIDs are not unique -> mistake in the selection process?")

# Split overlapping regions into 10 subsets ----
Nsplits <- 20
splits <- data.table(ov= unique(cmb$ov))
set.seed(1)
splits[, split:= sample(seq(Nsplits), .N, replace = T)]
cmb[splits, split:= i.split, on= "ov"]

# For each fold, define non-overlapping VISTA train/valid/test sets ----
folds <- sprintf("fold%02d", seq(Nsplits))
cmb[, (folds):= lapply(seq(folds), function(i) {
  fcase(ov %in% ov[seqnames=="chr18"], "test.shared", # Chr 18 always kept in test set
        split==i, "test",
        is.na(seqnames), "training", # 51 human-specific always kept in training
        split==(i+1 %% 20), "validation",
        default = "training")
})]

# Save ----
saveRDS(cmb, "db/folds/bulkATAC_folds.rds")