setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
require(AUCell)
require(Matrix)

# Import data ----
dat <- readRDS("db/single_cell/subsetted_sc_dataset.rds")

# Select cell clusters of interest ----
# dat$counts <- dat$counts[,dat$cells$cluster!="Other"]
# dat$cells <- dat$cells[cluster!="Other"]

# Compute cell rankings and auc ----
cells_rankings <- AUCell_buildRankings(dat$counts, plotStats = FALSE, verbose = TRUE)

# Import gene clusters ----
TFs <- readRDS("Rdata/motif_clusters_paper_3_tissues.rds")
TFs <- TFs[, .(mouse_id= unlist(mouse_id)), cluster]
TFs <- unique(na.omit(TFs[mouse_id %in% rownames(dat$counts)]))
TFs <- split(TFs$mouse_id, TFs$cluster)
names(TFs) <- paste0("Cluster ", names(TFs))
TFs <- c(TFs, list("Other TFs"= setdiff(rownames(dat$counts), unlist(TFs))))
auc <- AUCell_calcAUC(
  TFs,
  cells_rankings,
  aucMaxRank = nrow(dat$counts),
  nCores = getDTthreads()-1
)

# Retrieve AUCell scores ----
aucell_scores <- as.data.table(getAUC(auc), keep.rownames= "TF.cl")
aucell_scores <- melt(aucell_scores, id.vars= "TF.cl")
aucell_scores <- aucell_scores[TF.cl != "Other TFs"] # Only keep clustered TFs for now
aucell_scores[, cell.cl:= dat$cells$cluster[variable]]
aucell_scores[, zscore:= scale(value), TF.cl]

# Cast ----
aucell_scores[, TF.cl:= factor(TF.cl, c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"))]
aucell_scores[cell.cl=="Other", cell.cl:= "Other cell types"]
aucell_scores[, cell.cl:= factor(cell.cl, c("Cardiac muscle lineages", "Limb mesenchyme", "Neural tube", "Brain", "Other cell types"))]
mat <- dcast(aucell_scores, TF.cl~cell.cl, value.var = "zscore", fun.aggregate = mean)
mat <- as.matrix(mat, 1)

# Plot
vl_heatmap(mat,
           cluster.rows = F,
           legend.title = "AUCell z-score",
           pdf.file = "pdf/0_paper/heatmap_AUCell_zscore_motif_clusters_sc.pdf",
           breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3),
           pdf.cell.size = .24)
