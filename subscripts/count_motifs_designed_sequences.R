setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

before <- Biostrings::readDNAStringSet("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/evolutional_design/final_syth_enhancers/top30_tissue_20250220_step0.fa")
after <- Biostrings::readDNAStringSet("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/evolutional_design/final_syth_enhancers/top30_tissue_20250220.fa")
after <- after[match(names(before), names(after))]
# motifs <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/JASPAR_db/JASPAR2024_CORE_non-redundant_pwms_jaspar.rds")
motifs <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/non-redundant_pfms_jaspar_shenzhi.rds")

# Count motifs
if(!file.exists("db/counts/motifs/before.rds"))
{
  bCount <- vl_motif_counts(sequences = as.character(before),
                            pwm_log_odds = motifs,
                            genome = "mm10",
                            p.cutoff = 1e-4)
  saveRDS(bCount,
          "db/counts/motifs/before.rds")
}else
  bCount <- readRDS("db/counts/motifs/before.rds")

aCount <- if(!file.exists("db/counts/motifs/after.rds"))
{
  aCount <- vl_motif_counts(sequences = as.character(after),
                            pwm_log_odds = motifs,
                            genome = "mm10",
                            p.cutoff = 1e-4)
  saveRDS(aCount,
          "db/counts/motifs/after.rds")
}else
  aCount <- readRDS("db/counts/motifs/after.rds")

# Diff matrix
mat <- as.matrix(aCount-bCount)
rownames(mat) <- names(before)

# Tissues
tissues <- gsub("(.*)_.*", "\\1", names(before))

# Select columnes >= change
# sel <- mat
sel <- mat[, apply(mat, 2, function(x) sum(x!=0)>=3)]

# Cluster columns
col <- vl_heatmap(sel, cluster.rows = FALSE, plot = F)
clustered <- sel[, col$ccl$order]

# Plot
pdf("pdf/motif_counts_designed_sequences.pdf", width = 9, height = 8)
vl_par(mfrow= c(3,1),
       mai= c(0.1,0.1,0.5,1),
       omi= c(1,1,0,0))
for(tiss in unique(tissues)) {
  vl_heatmap(clustered[grepl(tiss, tissues), ],
             cluster.cols = FALSE,
             breaks = c(-5, 0, 5),
             main = tiss,
             show.colnames = tiss==tissues[length(tissues)],
             tilt.colnames = TRUE)
  
}
vl_par(mfrow= c(1,1))
vl_heatmap(sel,
           cluster.cols = FALSE,
           breaks = c(-5, 0, 5),
           main = tiss,
           show.colnames = tiss==tissues[length(tissues)],
           tilt.colnames = TRUE)
dev.off()