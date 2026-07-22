setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
mot <- readRDS("Rdata/annotated_PWMs.rds")[, .(motif, annot, cluster)]

# For each tissue ----
diff <- list()
for(tiss in c("heart", "limb", "midbrain")) {
  # Import counts and seq info
  counts <- readRDS(paste0("db/motifs/motif_counts_", tiss, "_sequences_0.0001_subject_bg.rds"))
  seq.info <- readRDS(paste0("Rdata/final_designed_enhancer_sequences_", tiss, ".rds"))
  .d <- counts[seq.info$label=="ledidi_12_14",]-counts[seq.info$label=="ledidi_12_14_ini",]
  .d[, selected:= seq.info[label=="ledidi_12_14", selected]]
  diff[[tiss]] <- .d
}
dat <- rbindlist(diff, idcol = "tissue")

# Melt ----
.m <- melt(dat, id.vars = c("tissue", "selected"), variable.name = "motif")
.m <- .m[, .(mean= mean(value)), .(tissue, motif, selected)]
.m <- merge(.m, mot, by= "motif")
sel <- .m[, motif[which.max(mean)], cluster]$V1
.m <- .m[motif %in% sel & !is.na(cluster)]

# Dcast ----
mat <- dcast(.m, annot+cluster~tissue+selected, value.var = "mean")
annot <- mat$annot
mat <- as.matrix(mat[, -1], 1)
sel <- apply(mat, 1, function(x) any(x>.5))
mat <- mat[(sel),]
colnames(mat) <- gsub("_FALSE", ": did not pass filtering", colnames(mat))
colnames(mat) <- gsub("_TRUE", " : passed filtering", colnames(mat))

# Plot heatmap ----
pdf("pdf/compare_designed_sequences_per_tissue.pdf", width = 3.75, height = 8)
vl_par(mai= c(1.3, 1.2, .4, 1.5))
br <- seq(-1, 4, .1)
col <- circlize::colorRamp2(c(min(br), 0, max(br)), colors = c("cornflowerblue", "white", "tomato"))(br)
vl_heatmap(mat,
           breaks = br,
           col = col,
           show.numbers = round(mat, 1),
           row.annotations = annot[sel],
           legend.title = "Mean counts",
           numbers.cex = .5,
           row.annotations.title = "Motif function")
dev.off()