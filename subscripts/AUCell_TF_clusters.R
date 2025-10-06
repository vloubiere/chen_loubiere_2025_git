setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
require(AUCell)
require(Matrix)

# Import data
dat <- readRDS("db/single_cell/subsetted_sc_dataset.rds")

# Compute cell rankings and auc
cells_rankings <- AUCell_buildRankings(dat$counts, plotStats = FALSE, verbose = TRUE)

# Import gene clusters
genes <- readRDS("Rdata/motif_clusters_paper_3_tissues.rds")
TF_clusters <- split(dat$genes$gene_name, dat$genes$cluster)
auc <- AUCell_calcAUC(TF_clusters, cells_rankings, aucMaxRank = nrow(dat$counts), nCores = getDTthreads()-1)

# Retrieve AUCell scores
aucell_scores <- as.data.table(getAUC(auc), keep.rownames= "TF.cl")
aucell_scores <- melt(aucell_scores, id.vars= "TF.cl")
aucell_scores[, cell.cl:= dat$cells$cluster[variable]]
aucell_scores[, zscore:= scale(value), TF.cl]

# Cast
mat <- dcast(aucell_scores, TF.cl~cell.cl, value.var = "zscore", fun.aggregate = mean)
mat <- as.matrix(mat, 1)
rownames(mat) <- paste0("TF cluster ", rownames(mat))
rownames(mat)[rownames(mat)=="TF cluster 0"] <- "Other TFs"
colnames(mat)[colnames(mat)=="Other"] <- "Other cell types"
mat <- mat[, c("Brain", "Neural tube", "Limb mesenchyme", "Cardiac muscle lineages", "Other cell types")]

# Plot
vl_heatmap(mat,
           legend.title = "AUCell z-score",
           pdf.file = "pdf/0_paper/heatmap_AUCell_zscore_motif_clusters_sc.pdf",
           breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
           pdf.cell.size = .24)
