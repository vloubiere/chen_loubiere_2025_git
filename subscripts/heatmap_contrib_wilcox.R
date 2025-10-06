setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import data ----
dat <- readRDS("db/contributions/mean_contrib_per_motif_instance.rds")

FC <- dat[, .(
  FC= mean(contrib.mot)/mean(contrib.ctl, na.rm= T),
  p.value= wilcox.test(contrib.mot, contrib.ctl)$p.value
), .(tissue, dataset, motif)]
FC[, c("cluster", "motif"):= tstrsplit(motif, "__")]
res <- FC[, .(mean.FC= mean(FC), pval= metap::sumlog(p.value)$p, motif= motif[which.min(p.value)]), .(dataset, cluster, tissue)]
res[, padj:= p.adjust(pval, "fdr")]

# Matrix ----
mat <- dcast(res, cluster~tissue+dataset, value.var = "mean.FC")
mat <- as.matrix(mat, 1)


vl_heatmap(mat, breaks= seq(-10, 10, length.out= 21))
vl_heatmap(scale(mat), breaks= seq(-3, 3, length.out= 21))

# Subset and dcast ----
sel <- dat[, .SD[which.max(zscore)], cluster][zscore>=zsore_cutoff]$motif
scaled <- dcast(dat[sel, on= "motif"], id~tissue+dataset, value.var = "zscore")
scaled <- as.matrix(scaled, 1)

# Clip outliers
clipped <- apply(scaled, 2, function(x) {
  lim <- quantile(x, c(0.05, 0.95))
  x[x<lim[1]] <- lim[1]
  x[x>lim[2]] <- lim[2]
  x
})

# Plot
col.cl <- c("Heart", "Heart", "Limb", "Limb", "Midbrain", "Midbrain")
col.cl <- factor(col.cl, c("Midbrain", "Heart", "Limb"))
pdf("pdf/0_paper/heatmap_scaled_motif_contrib_3_tissues.pdf", height = 8, width = 7)
vl_par(mai= c(1,3,.5,1.5))
hm <- vl_heatmap(clipped,
                 cutree.rows = 4,
                 cluster.cols = col.cl,
                 breaks = seq(-3.5, 3.5, length.out= 21),
                 legend.title = "z-score",
                 row.gap.width = .5)
dev.off()

hm <- hm$rows
hm[dat, TFs:= i.TFs, on= "name==id"]
saveRDS(hm, "Rdata/motif_clusters_paper_3_tissues.rds")
