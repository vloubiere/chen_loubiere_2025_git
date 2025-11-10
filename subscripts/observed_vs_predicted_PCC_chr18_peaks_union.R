setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import peaks ----
files <- c("db/peaks/ATAC/midbrain_peaks.narrowPeak",
           "db/peaks/ATAC/limb_peaks.narrowPeak",
           "db/peaks/ATAC/heart_peaks.narrowPeak")
peaks <- lapply(files, importBed)
peaks <- rbindlist(peaks)
peaks <- collapseBed(peaks[seqnames=="chr18"])

# Compute signal ----
tracks <- c("db/bw/observed/midbrain_treat_pileup.bigwig",
            "db/bw/observed/heart_treat_pileup.bigwig",
            "db/bw/observed/limb_treat_pileup.bigwig",
            "db/bw/predicted/midbrain_predicted_accessibility_chr18.bw",
            "db/bw/predicted/heart_predicted_accessibility_chr18.bw",
            "db/bw/predicted/limb_predicted_accessibility_chr18.bw")
sig <- lapply(tracks, bwCoverage, bed= peaks)
mat <- do.call(cbind, sig)
colnames(mat) <- c("Midbrain.obs", "Heart.obs", "Limb.obs",
                   "Midbrain.pred", "Heart.pred", "Limb.pred")

# Compute correlations ----
cor <- data.table(
  tissue= c("Midbrain", "Heart", "Limb"),
  PCC= c(
    cor(log2(mat[, "Midbrain.pred"]+1), log2(mat[, "Midbrain.obs"]+1)),
    cor(log2(mat[, "Heart.pred"]+1), log2(mat[, "Heart.obs"]+1)),
    cor(log2(mat[, "Limb.pred"]+1), log2(mat[, "Limb.obs"]+1))
  )
)

# Save ----
saveRDS(cor,
        "db/PCC/obs_vs_pred_per_tissue_chr18_peaks_union.rds")

