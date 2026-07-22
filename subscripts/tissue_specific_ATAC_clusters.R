setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import ATAC peaks ----
ATAC <- readRDS("db/peaks/bulkENCODE_confident_ATAC_peaks.rds")[, 1:5]

# Import VISTA mm10 peaks ----
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
vista$class <- vista$genome <- vista$coor_hg38 <- NULL
vista <- vista[!is.na(coor_mm10)]
vista[, c("seqnames", "start", "end"):= importBed(coor_mm10)[, 1:3, with= F]]
vista$coor_mm10 <- NULL

# Restrict to overlapping regions ----
vista <- intersectBed(vista, ATAC)
ATAC <- intersectBed(ATAC, vista)
vista$strand <- NULL
ATAC$strand <- NULL

# Compute overlaps for each ATAC class (including 'all') ----
mVista <- melt(vista, id.vars = c("peakID", "seqnames", "start", "end"), value.name = "active")
mATAC <- rbind(
  data.table::copy(ATAC)[, class:= "all"],
  ATAC
)
ov <- overlapBed(mVista, mATAC)
ov[, tissue:= mVista$variable[idx.a]]
ov[, active:= mVista$active[idx.a]]
ov[, class:= mATAC$class[idx.b]]

# Plot ----
perc <- dcast(ov, tissue~class, value.var = "active", fun.aggregate = function(x) sum(x)/length(x)*100)
perc <- as.matrix(perc, 1)
numbers <- dcast(ov, tissue~class, value.var = "active", fun.aggregate = function(x) paste0("n=", sum(x)))
numbers <- as.matrix(numbers, 1)
numbers <- numbers[, colnames(perc)]
vl_heatmap(
  perc,
  show.numbers = round(perc, 1),
  cluster.rows = F,
  legend.title = "TPR",
  main= "Clustered peaks",
  pdf.file = "pdf/1_tests_review/overlap_VISTA_clustered_ATAC_peaks.pdf",
  pdf.cell.size = .3,
  pdf.close = FALSE
)
text(rep(seq(ncol(numbers)), each= nrow(numbers)), rep(seq(nrow(numbers)), ncol(numbers)), numbers, cex= .4, pos= 1, offset= 0.4)
title(xlab= "ATAC peaks", line = 5)
title(ylab= "VISTA labels", line = 3)
dev.off()