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
enr <- enr[motif %in% enr[rank <= 25, motif]]

# As matrix ----
log2OR <- dcast(enr, name~tissue+class, value.var = "log2OR")
log2OR <- as.matrix(log2OR, 1)
padj <- dcast(enr, name~tissue+class, value.var = "padj")
padj <- as.matrix(padj, 1)
log2OR[padj>0.05] <- 0
lim <- quantile(log2OR, c(0.05, 0.95))
clip <- log2OR
clip[clip<lim[1]] <- lim[1]
clip[clip>lim[2]] <- lim[2]

# Import motif counts at validate seq ----
mot <- readRDS("db/motifs/revision_motif_counts_validated_synthetic_enhancer_seq.rds")
mot.mat <- dcast(mot[variable %in% enr$name], variable~tissue+id, value.var = "value")
mot.mat <- as.matrix(mot.mat, 1)
mot.mat <- mot.mat[rownames(log2OR),]

# Simplify names ----
simp.names <- as.data.table(tstrsplit(rownames(log2OR), "_", keep= c(1,3)))[, paste0(V1, " (", V2, ")")]
simp.names <- gsub(".mouse", "", simp.names)
rownames(log2OR) <- simp.names
rownames(mot.mat) <- simp.names

# Motif counts validated enhancers ----
Nclust <- 6
pdf("pdf/_revision/heatmap_top_motif_enrich.pdf", width = 7.75, height = 8)
layout(matrix(c(1,2), ncol= 2), widths = c(9, 12))
vl_par(mai= c(.2, .2, .2, .2), omi= c(.9, 2.2, .9, 1.2), lwd= .5)
hm <- vl_heatmap(
  log2OR,
  cluster.rows = clip,
  # clustering.method = "ward.D",
  breaks = seq(-3, 3, length.out= 21),
  cluster.cols = T,
  cutree.rows = Nclust,
  show.row.clusters = F,
  show.legend = "top",
  legend.title = "OR (log2) vs. closed genomic seq."
)
vl_heatmap(
  mot.mat,
  cluster.rows = clip,
  breaks = seq(0, 8),
  cutree.rows = Nclust,
  show.row.clusters = F,
  show.rownames= F,
  cluster.cols = T,
  show.numbers = mot.mat,
  numbers.cex = .5,
  legend.title = "Motif count"
)
dev.off()
