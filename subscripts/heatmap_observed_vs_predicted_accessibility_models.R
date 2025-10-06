setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v2.rds")
meta <- meta[tissue %in% c("midbrain", "heart", "limb")]
meta <- meta[dataset=="accessibility" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set=="test"]

# Import predicted values ----
dat <- meta[, {
  # Import
  .c <- .SD[, {
    merge(
      fread(obs_file)[grepl("\\+_0__.*", ID)], # Keep center tiles only
      fread(pred_file),
      by.x= "ID",
      by.y= "location"
    )
  }, .(obs_file, pred_file, replicate, fold)]
  # Mean across fold/reps
  .c[, lapply(.SD, mean), ID, .SDcols= c("score", "Predictions")]
}, tissue]

# Split ID and tissue-specific label ----
dat[, c("ID", "label"):= tstrsplit(ID, "__")]
dat[, ID:= tstrsplit(ID, "_", keep= 1)]

# Compute distance closest promoter ----
prom <- rtracklayer::import("/groups/stark/vloubiere/projects/ORFTRAP_1/db/gtf/gencode.vM25.basic.annotation.gtf.gz")
prom <- as.data.table(prom)
tss <- resizeBed(prom[type=="gene", .(seqnames, start, end)], "start", 0, 0)
uniq.tiles <- data.table(ID= unique(dat$ID))
uniq.tiles[, c("seqnames", "start", "end"):= importBed(ID)[, .(seqnames, start, end)]]
dist <- closestBed(uniq.tiles, tss)
uniq.tiles$dist <- dist[, mean(abs(dist)), idx.a]$V1

# Observed Score matrix ----
score <- dcast(dat, ID~tissue, value.var = "score")
score <- as.matrix(score, 1)
# Clip outliers
score <- apply(score, 2, function(x) {
  clip <- quantile(x, .5, .95)
  x[x<clip[1]] <- clip[1]
  x[x>clip[2]] <- clip[2]
  x
})
# Scale (z-score)
scaled.score <- scale(score)

# Predicted matrix ----
pred <- dcast(dat, ID~tissue, value.var = "Predictions")
pred <- as.matrix(pred, 1)
# Scale (z-score)
scaled.pred <- scale(pred)

# Select only tiles that are active in at least one tissue ----
act.sel <- apply(scaled.score, 1, function(x) any(x>1))

# Plot ----
breaks <- seq(.5, 4, length.out= 21)
col <- colorRampPalette(c("white", "red"))(21)
Nclust <- 6

# Cluster ----
kcl <- vlite::vl_heatmap(
  scaled.score[act.sel,],
  kmeans.k = Nclust,
  plot = F
)$rows

# Plot
pdf("pdf/0_paper/heatmap_observed_vs_predicted_accessibility_models.pdf",
    width = 7,
    height = 4)
layout(matrix(1:3, nrow= 1), widths = c(.8, .8, 1.8))
vl_par(mai= c(.2, .2, .2, .2),
       omi= c(1, 1, 1, 1))
hm <- vlite::vl_heatmap(
  scaled.score[act.sel,],
  cluster.rows = factor(kcl$cluster, c("4", "2", "6", "3", "5", "1")),
  breaks = breaks,
  col= col,
  show.rownames = F,
  cluster.cols= F,
  show.legend = "top",
  legend.title = "ATAC-seq z-score",
  legend.cex = .6,
  show.row.clusters = "left"
)
vlite::vl_heatmap(
  scaled.pred[act.sel, hm$cols[order(x.pos), name]],
  cluster.rows = factor(kcl$cluster, c("4", "2", "6", "3", "5", "1")),
  cluster.cols= F,
  kmeans.k = Nclust,
  breaks = breaks,
  col= col,
  show.rownames = F,
  show.legend = "top",
  legend.cex = .6,
  legend.title = "Predicted z-score",
  show.row.clusters = F
)
# Add distance to TSS per cluster 
cl.dist <- merge(hm$rows, uniq.tiles, by.x= "name", by.y= "ID")
vl_par(mai= c(.5, 0.5, .5, .2),
       mgp= c(1.5, .35, 0),
       cex.axis= .4)
vl_boxplot(
  dist~cluster,
  cl.dist,
  xlab= "Clusters",
  ylab= "Distance to closest TSS"
)
dev.off()
