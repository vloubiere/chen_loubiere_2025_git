setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# Import mean contrib per motif per rep/fold
dat <- readRDS("/groups/stark/shenzhi.chen/projects/mouse_enhancer_paper/1st_revision/Rdata/all_motif_contri_score.rds")

# Z-score
dat[, zscore:= scale(contrib), .(dataset, tissue, replicate, fold)]

# Add heatmap motifs information
hm <- readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/motif_clusters_paper_3_tissues.rds")
setorderv(hm, "cluster")
dat <- dat[id %in% hm$name]
dat[hm, cluster:= i.cluster, on= "id==name"]
dat[, id:= factor(id, unique(hm$name))]
setorderv(dat, "id")
dat[, motif:= factor(motif, unique(motif))]

# Merge accessibility and activity values
dat <- merge(
  dat[dataset=="accessibility", .(motif, cluster, tissue, replicate, fold, zscore)],
  dat[dataset=="activity", .(motif, tissue, replicate, fold, zscore)],
  by= c("motif", "tissue", "replicate", "fold"),
  suffixes= c(".acc", ".act")
)
dat[, at:= .GRP, motif]

# Plot ----
width <- .4
Cc <- c("white", "lightgrey")
pdf("pdf/_revision/boxplot_contrib_TL.pdf", width = 8, height = 6.5)
vl_par(mai= c(.9, .9, .4, .9), mfrow= c(3,1))
dat[, {
  Nmot <- length(unique(motif))
  # Barplot accessibility
  vl_barplot(
    split(zscore.acc, motif),
    pch.jitter= NULL,
    at= seq(length(unique(motif)))-width/2,
    width = width,
    col= Cc[1],
    pch.col= "black",
    ylim= range(c(zscore.acc, zscore.act)),
    xaxt= "n",
    ylab= "Contribution (z-score)",
    main= tissue
  )
  title(xlab= "Motifs", line = 5)
  # Barplot activity
  vl_barplot(
    split(zscore.act, motif),
    pch.jitter= NULL,
    at= seq(length(unique(motif)))+width/2,
    width = width,
    col= Cc[2],
    add= T,
    pch.col= "black",
    xaxt= "n"
  )
  # Axis
  tiltAxis(
    x = seq(length(unique(motif))),
    labels = levels(motif)
  )
  # Connect points
  .SD[, {
    segments(
      at[1]-width/2,
      zscore.acc[1],
      at[1]+width/2,
      zscore.act[1],
      lwd= .5,
      col= "grey20"
    )
  }, .(motif, at, zscore.acc, zscore.act)]
  # Add motif clusters
  .SD[, {
    segments(
      min(at)-width/2,
      par("usr")[4],
      max(at)+width/2,
      par("usr")[4]
    )
    text(
      mean(range(at)),
      par("usr")[4],
      paste0("Cluster ", cluster),
      pos= 3,
      xpd= T,
      offset= 0.2
    ) 
  }, cluster]
  vl_legend(
    fill= Cc,
    legend= c("Accessibility", "Activity"),
    border = "black"
  )
  print("")
}, tissue]
dev.off()