setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import enrichment ----
enr <- readRDS("db/motifs/revision_motif_enrich_atac_vista_designed_vs_rdm_genomic.rds")
enr[, sig:= padj<0.05 & log2OR>0 & set_hit>=5]
setorderv(enr, c("sig", "log2OR"), -1)

# Select top motif per cluster ----
enr[, cluster:= tstrsplit(name, "__", keep= 1)]
enr <- enr[, .SD[motif==motif[1]], cluster]

# Select top motifs per class ----
enr[(sig), rank:= rowid(class, tissue)]
enr <- enr[motif %in% enr[rank <= Inf, motif]]

# As matrix ----
log2OR <- dcast(enr, cluster~tissue+class, value.var = "log2OR")
log2OR <- as.matrix(log2OR, 1)
padj <- dcast(enr, cluster~tissue+class, value.var = "padj")
padj <- as.matrix(padj, 1)
log2OR[padj>0.05] <- 0
lim <- quantile(log2OR, c(0.05, 0.95))
clip <- log2OR
clip[clip<lim[1]] <- lim[1]
clip[clip>lim[2]] <- lim[2]

# Import motif counts at validate seq ----
mot <- readRDS("db/motifs/revision_motif_counts_validated_synthetic_enhancer_seq.rds")
mot.mat <- dcast(mot[variable %in% enr$name], variable~tissue+id, value.var = "value")
mot.mat[, variable:= tstrsplit(variable, "__", keep= 1)]
mot.mat <- as.matrix(mot.mat, 1)
mot.mat <- mot.mat[rownames(log2OR),]

# Motif counts validated enhancers ----
Nclust <- 7
pdf("pdf/_revision/heatmap_all_motifs_enrich.pdf", width = 3.2, height = 4)
vl_par(lwd= .5)
hm <- vl_heatmap(
  log2OR,
  cluster.rows = clip,
  breaks = seq(-3, 3, length.out= 21),
  cluster.cols = T,
  cutree.rows = Nclust,
  show.row.clusters = F,
  show.rownames = F,
  show.legend = "top",
  legend.title = "OR (log2) vs. closed genomic seq."
)
dev.off()
