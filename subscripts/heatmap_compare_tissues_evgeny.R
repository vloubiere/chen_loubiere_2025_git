setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
mot <- readRDS("Rdata/annotated_PWMs.rds")[, .(motif, annot, cluster)]

# For each tissue ----
if(!exists("dat") || !is.data.table(dat)) {
  dat <- list()
  for(tiss in c("heart", "limb", "midbrain")) {
    # Import counts and seq info
    counts <- readRDS(paste0("db/motifs/motif_counts_", tiss, "_sequences_0.0001_subject_bg.rds"))
    seq.info <- readRDS(paste0("Rdata/final_designed_enhancer_sequences_", tiss, ".rds"))
    .d <- counts[seq.info$label=="ledidi_12_14",]-counts[seq.info$label=="ledidi_12_14_ini",]
    seq.info <- seq.info[label=="ledidi_12_14", .(blast, active, specific, selected)]
    .d <- cbind(seq.info, .d)
    dat[[tiss]] <- .d
  }
  dat <- rbindlist(dat, idcol = "tissue")
}

# Melt ----
.m <- melt(dat,
           id.vars = c("tissue", "blast", "active", "specific", "selected"),
           variable.name = "motif")
.m[, All:= TRUE]
.m[, noBlast:= !blast]
.m$blast <- NULL

# Compute means ----
counts <- melt(.m, id.vars = c("tissue", "motif", "value"), value.name = "sel")
counts <- counts[(sel)]
counts <- counts[, .(mean= mean(value)), .(tissue, motif, variable)]

# Select motif with higher counts per cluster  ----
counts[mot, c("cluster", "annot"):= .(i.cluster, i.annot), on= "motif"]
sel <- counts[, motif[which.max(mean)], cluster]$V1
sel <- counts[motif %in% sel & !is.na(cluster)]
sel[, variable:= factor(variable, c("All", "selected", "noBlast", "active", "specific"))]

# Dcast ----
mat <- dcast(sel, annot+cluster~tissue+variable, value.var = "mean")
annot <- mat$annot
mat <- as.matrix(mat[, -1], 1)
sel <- apply(mat, 1, function(x) any(round(x, 1)>=.5))
mat <- mat[(sel),]

# Diff heatmap ----
diff.hm <- cbind(
  mat[,"heart_selected", drop= F]-mat[, "heart_All", drop= FALSE],
  mat[,"limb_selected", drop= F]-mat[, "limb_All", drop= FALSE],
  mat[,"midbrain_selected", drop= F]-mat[, "midbrain_All", drop= FALSE]
  )
colnames(diff.hm) <- unlist(tstrsplit(colnames(diff.hm), "_", keep= 1))

# Plot heatmap ----
pdf("pdf/compare_designed_sequences_per_tissue_evegeny.pdf", width = 12, height = 9.5)
vl_par(mai= c(1.3, 1, .4, 2.5),
       mfrow= c(1,2))
br <- seq(-1, 4, .1)
Cc <- c("cornflowerblue", "white", "tomato")
col <- circlize::colorRamp2(c(min(br), 0, max(br)), colors = Cc)(br)
vl_heatmap(mat,
           breaks = br,
           col = col,
           show.numbers = round(mat, 1),
           row.annotations = annot[sel],
           cluster.cols = unlist(tstrsplit(colnames(mat), "_", keep= 1)),
           legend.title = "Mean counts (ledidi-init.seq)",
           numbers.cex = .5,
           row.annotations.title = "Manual annotation")
par(mai= c(1.3, 1, .4, 4.5))
br <- seq(-.5, 2, .1)
col <- circlize::colorRamp2(c(min(br), 0, max(br)), colors = Cc)(br)
vl_heatmap(diff.hm,
           breaks = br,
           cluster.rows = mat,
           col = col,
           show.numbers = round(diff.hm, 1),
           row.annotations = annot[sel],
           legend.title = "Difference (Selected-All)",
           numbers.cex = .5,
           row.annotations.title = "Manual annotation")
dev.off()