setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
require(Seurat)
require(Matrix)

# Import data
mat <- fread(
  "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/MouseAtlas/GSE119945_gene_count.txt.gz",
  sel= 1:3,
  col.names = c("gene", "cell", "count")
)

# Subset E11.5 cells
cells <- fread("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/MouseAtlas/GSE119945_cell_annotate.csv.gz")
sub <- mat[cell %in% which(cells$day=="11.5")]
sub[, cell:= .GRP, cell]
cells <- cells[day=="11.5"]

# Retrieve gene names
genes <- fread("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/MouseAtlas/GSE119945_gene_annotate.csv.gz")
genes[, gene_id:= tstrsplit(gene_id, "[.]", keep= 1)]
sub[, gene_name:= genes$gene_short_name[gene]]
sub[, gene_id:= genes$gene_id[gene]]

# Subset TF genes
all_TFs <- readRDS("/groups/stark/vloubiere/motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")
clean <- sub[gene_id %in% na.omit(unlist(all_TFs$meta$mouse_id))]

# Subset TF genes
clean[, gene_id:= factor(gene_id)]
s.mat <- Matrix::sparseMatrix(i= clean$gene_id, j= clean$cell, x = clean$count)
rownames(s.mat) <- levels(clean$gene_id)
colnames(s.mat) <- seq(max(clean$cell))

# Retrieve cell clusters
cl <- fread("db/single_cell/clusters.txt")
cells[, cluster:= cl$simp_cluster[Main_Cluster]]
cells[, total.cl:= .N, cluster]
cells.cl <- cells[, .(cluster, total.cl, tsne_1, tsne_2)]

# Remove missing cells
empty <- which(colSums(s.mat)==0)
s.mat <- s.mat[, -c(empty)]
cells.cl <- cells.cl[-c(empty)]

# Remove NA cells (no clusters)
s.mat <- s.mat[, !is.na(cells.cl$cluster)]
cells.cl <- cells.cl[!is.na(cluster)]

# Save ----
saveRDS(list(counts= s.mat, cells= cells.cl), "db/single_cell/subsetted_sc_dataset.rds")
