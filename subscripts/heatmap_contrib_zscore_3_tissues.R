setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import data ----
dat <- readRDS("db/contributions/mean_contrib_per_motif_instance.rds")
dat <- dat[tissue %in% c("midbrain", "limb", "heart"), .(contrib= mean(contrib.mot)), .(motif, dataset, tissue)]
dat[, id:= motif]
dat[, c("cluster", "motif"):= tstrsplit(motif, "__")]

# Scale and select high contrib motifs ----
dat[, zscore:= scale(contrib), .(dataset, tissue)]
zscore_cutoff <- 2
sel <- dat[id %in% dat[zscore>=zscore_cutoff, id]]

# Get motif names and retain only best motif per cluster----
all_TFs <- readRDS("/groups/stark/vloubiere/motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")
sel[, TFs:= , cluster]
sel[, TFs:= paste0(sort(unique(unlist(tstrsplit(unique(TFs), "\\+")))), collapse= ","), cluster]

# Subset and dcast ----
sel.mot <- sel[, .SD[which.max(zscore)], cluster][zscore >= zscore_cutoff]$motif
scaled <- dcast(sel[sel.mot, on= "motif"], id~tissue+dataset, value.var = "zscore")
scaled <- as.matrix(scaled, 1)

# Clip outliers ----
clipped <- apply(scaled, 2, function(x) {
  lim <- quantile(x, c(0.05, 0.95))
  x[x<lim[1]] <- lim[1]
  x[x>lim[2]] <- lim[2]
  x
})

# clustering ----
set.seed(3453)
hc <- hclust(dist(clipped))
cl <- cutree(hc, 4)
cl <- c(1, 4, 2, 3)[cl]

# Plot ----
pdf("pdf/0_paper/heatmap_scaled_motif_contrib_3_tissues.pdf", height = 8, width = 7)
vl_par(mai= c(1,3,.5,1.5))
hm <- vl_heatmap(clipped,
                 cluster.rows = cl,
                 cluster.seed = 4,
                 cluster.cols = gsub("_accessibility$|_activity$", "", colnames(clipped)),
                 breaks = seq(-3.5, 3.5, length.out= 21),
                 legend.title = "z-score",
                 row.gap.width = .5)
dev.off()

hm <- hm$rows
hm[sel, TFs:= i.TFs, on= "name==id"]
saveRDS(hm, "Rdata/motif_clusters_paper_3_tissues.rds")
