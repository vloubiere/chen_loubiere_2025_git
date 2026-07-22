setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import ATAC peaks ----
peak.files <- list.files("db/peaks/ATAC/", full.names = T)
peaks <- lapply(peak.files, importBed)
names(peaks) <- unlist(tstrsplit(basename(peak.files), "_", keep= 1))
peaks <- rbindlist(peaks, idcol = "tissue")
peaks <- collapseBed(peaks[signalValue>3 & qValue>5])

# Import RNA coordinates ----
RNA <- readRDS("Rdata/tissue_specific_genes_RNA.rds")

# Assign peaks to genes ----
cl <- closestBed(peaks, RNA)
plot(density(abs(cl$dist)), xlim= c(-100000, 500000))
abline(v= 200000)
cl <- cl[abs(dist)<200000]

# Define tissue-specific peaks ----
peaks.cl <- peaks[cl$idx.a]
peaks.cl[, tissue:= RNA[cl$idx.b, class]]

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
cols <- c("heart-specific",
          "limb-specific",
          "facialProminence-specific", 
          "forebrain-specific",
          "midbrain-specific",
          "hindbrain-specific",
          "neuralTube-specific",
          "liver-specific",
          "heart-enriched", 
          "limb-enriched",
          "forebrain-enriched",
          "midbrain-enriched",
          "hindbrain-enriched",
          "neuralTube-enriched",
          "liver-enriched",
          "Shared")
stopifnot(all(cols %in% colnames(perc)) & all(colnames(perc) %in% cols))
perc <- perc[, cols]

# Compute Nominator ----
numbers <- dcast(ov, tissue~class, value.var = "active", fun.aggregate = function(x) paste0("n=", sum(x)))
numbers <- as.matrix(numbers, 1)
numbers <- numbers[, colnames(perc)]

# Plot ----
vl_heatmap(
  perc,
  show.numbers = round(perc, 1),
  cluster.rows = F,
  legend.title = "TPR",
  main= "Overlapped peaks",
  pdf.file = "pdf/1_tests_review/overlap_VISTA_specific_closest_gene.pdf",
  pdf.cell.size = .3,
  pdf.close = FALSE
)
text(rep(seq(ncol(numbers)), each= nrow(numbers)), rep(seq(nrow(numbers)), ncol(numbers)), numbers, cex= .4, pos= 1, offset= 0.4)
title(xlab= "ATAC peaks", line = 3)
title(ylab= "VISTA labels", line = 3)
dev.off()
