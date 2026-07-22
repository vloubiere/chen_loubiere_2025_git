setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import K27Ac/K4me1 peaks ----
folder <- "/groups/stark/shenzhi.chen/projects/accessibility_model_enhancer_design_17112025/"
meta <- readRDS(paste0(folder, "Rdata/mouse_e11.5_ENCODE_20251204/metadata.rds"))
meta[, files:= paste0(folder, files)]
meta <- meta[dataset %in% c("H3K27ac", "H3K4me1")]
dat <- meta[, {
  .c <- fread(files, sel= c(1,2,3,7,9), col.names = c("seqnames", "start", "end", "signalValue", "qValue"))
  .c[, start:= start+1]
}, .(dataset, tissue)]

# Only retain K27Ac peaks that have K4me1 ----
peaks <- dat[, {
  intersectBed(.SD[dataset=="H3K27ac"], .SD[dataset=="H3K4me1"])
}, tissue]
peaks <- peaks[signalValue>2 & qValue>2]

# Define tissue-specific peaks ----
ov <- overlapBed(peaks, peaks)
ov[, tissue.a:= peaks$tissue[idx.a]]
ov[, tissue.b:= peaks$tissue[idx.b]]
ov.idx <- unique(ov[tissue.a!=tissue.b, idx.a])
peaks[, tissueSpecific:= ifelse(.I %in% ov.idx, FALSE, TRUE)]
peaks.cl <- rbind(
  peaks[(tissueSpecific)],
  collapseBed(peaks[(!tissueSpecific)])[, tissue:= "shared"],
  collapseBed(peaks)[, tissue:= "allPeaks"],
  fill= T
)
peaks.cl <- peaks.cl[, .(seqnames, start, end, tissue)]
peaks.cl[, .N, tissue]

# Import VISTA mm10 peaks ----
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
vista$class <- vista$genome <- vista$coor_hg38 <- NULL
vista <- vista[!is.na(coor_mm10)]
vista[, c("seqnames", "start", "end"):= importBed(coor_mm10)[, 1:3, with= F]]
vista$coor_mm10 <- NULL

# Restrict to overlapping regions ----
vista <- intersectBed(vista, peaks.cl)
peaks.cl <- intersectBed(peaks.cl, vista)
vista$strand <- NULL
peaks$strand <- NULL

# Compute overlaps for each ATAC class (including 'all') ----
mVista <- melt(vista, id.vars = c("peakID", "seqnames", "start", "end"), value.name = "active")
ov <- overlapBed(mVista, peaks.cl)
ov[, tissue:= mVista$variable[idx.a]]
ov[, active:= mVista$active[idx.a]]
ov[, class:= peaks.cl$tissue[idx.b]]

# Compute percentage and order columns ----
perc <- dcast(ov, tissue~class, value.var = "active", fun.aggregate = function(x) sum(x)/length(x)*100)
perc <- as.matrix(perc, 1)
cols <- c("heart", "limb", "forebrain", "midbrain", "hindbrain", "neuralTube", "liver", "shared", "allPeaks")
stopifnot(all(cols %in% colnames(perc)) & all(colnames(perc) %in% cols))
perc <- perc[, cols]

# Compute Nominator ----
numbers <- dcast(ov, tissue~class, value.var = "active", fun.aggregate = function(x) paste0("n=", sum(x)))
numbers <- as.matrix(numbers, 1)
numbers <- numbers[, colnames(perc)]

# Plot
vl_heatmap(
  perc,
  show.numbers = round(perc, 1),
  cluster.rows = F,
  legend.title = "TPR",
  main= "Overlapped peaks",
  pdf.file = "pdf/1_tests_review/overlap_VISTA_K27Ac_K4me1_peaks.pdf",
  pdf.cell.size = .3,
  pdf.close = FALSE
)
text(rep(seq(ncol(numbers)), each= nrow(numbers)), rep(seq(nrow(numbers)), ncol(numbers)), numbers, cex= .4, pos= 1, offset= 0.4)
title(xlab= "ATAC peaks", line = 2)
title(ylab= "VISTA labels", line = 3)
dev.off()