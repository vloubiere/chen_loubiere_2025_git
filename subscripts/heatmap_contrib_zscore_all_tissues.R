setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import data ----
dat <- readRDS("db/contributions/mean_contrib_per_motif_instance.rds")
dat <- dat[, .(contrib= mean(contrib.mot)), .(motif, dataset, tissue)]
dat[, id:= motif]
dat[, c("cluster", "motif"):= tstrsplit(motif, "__")]

zsore_cutoff <- 2

# Scale and dcast ----
dat[, zscore:= scale(contrib), .(dataset, tissue)]
dat[zscore >= zsore_cutoff, TFs:= tstrsplit(id, "_", keep= 3), cluster]
dat[, TFs:= paste0(sort(unique(na.omit(TFs))), collapse = ","), cluster]

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
col.cl <- rep(c("Forebrain", "Heart", "Hindbrain", "Limb", "Midbrain", "Neural tube"), each= 2)
col.cl <- factor(col.cl, c("Midbrain", "Neural tube", "Forebrain", "Hindbrain", "Heart", "Limb"))

pdf("pdf/0_paper/heatmap_scaled_motif_contrib_all_tissues.pdf", height = 8, width = 8)
vl_par(mai= c(1,3,.5,1.5))
hm <- vl_heatmap(clipped,
                 cutree.rows = 6,
                 cluster.cols= col.cl,
                 breaks = seq(-3.5, 3.5, length.out= 21),
                 legend.title = "z-score",
                 row.gap.width = .5)
dev.off()

hm <- hm$rows
hm[dat, TFs:= i.TFs, on= "name==id"]
saveRDS(hm, "Rdata/motif_clusters_paper_all_tissues.rds")
