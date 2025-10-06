setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import data ----
dat <- readRDS("db/contributions/mean_contrib_per_motif_instance.rds")
dat <- dat[tissue %in% c("midbrain", "limb", "heart")]
dat <- dat[, .(
  contrib= mean(contrib.mot), 
  pval= wilcox.test(contrib.mot, contrib.ctl, alternative="greater")$p.value
), .(motif, dataset, tissue)]
dat[, padj:= p.adjust(pval, "fdr"), .(motif, dataset, tissue)]
dat[, id:= motif]
dat[, c("cluster", "motif"):= tstrsplit(motif, "__")]

# Get TFs ----
TFs <- readRDS("/groups/stark/vloubiere/motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")$meta
dat[TFs, TFs:= i.TF, on= "motif==Motif"]
dat[, zscore:= scale(contrib), .(dataset, tissue)]
zscore_cutoff <- 2
dat[zscore >= zscore_cutoff & padj<0.001, TFs:= .(.(sort(unique(unlist(TFs))))), cluster]

# Subset and dcast ----
sel <- dat[, .SD[which.max(zscore)], cluster][zscore >= zscore_cutoff]$motif
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
pdf("pdf/0_paper/heatmap_scaled_motif_contrib_3_tissues_single_TF.pdf", height = 8, width = 7)
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
saveRDS(hm, "Rdata/motif_clusters_paper_3_tissues_single_TF.rds")
