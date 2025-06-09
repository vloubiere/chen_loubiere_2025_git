setwd("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Input ----
tissue <- "heart"
enr.file <- "Rdata/enhancer_design/evolutional_design/motif_enrichment/vista_enhancer_with_all_motif_enrich.rds"
contrib.file <- "result/model_evaluation/motifs/jeff_relative_weight_matrix_enhancer_model1_bulkATAC_10xAug_nogoBal_noW.rds"
# mot.counts.file <- "db/counts/motifs/heart_all_sequence_counts.rds"
mot.counts.file <- "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/motifs/motif_counts_genome_bg.rds"
# mot.counts.file <- "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/motifs/motif_counts_subject_bg.rds"
# Row labels
seq.info <- "Rdata/motifs_enrichment_analysis/all_sequences_for_motifs_analysisheart.rds"
label.column <- "label"
label.lvls <- c("rdm", "dinuc.adj", "vista_inactive_all_tissues",
                "designed_heart", "ledidi", "vista_ts",  "vista_other_tissues")
label_colors <- c("grey90", "grey60", "black", "limegreen", "lightgreen", "red", "blue")
contrib.cutoffs <- c(0.0002, 0.0004)
padj.cutoffs <- c(5e-2, 1e-7)

# Retrieve motif clusters ----
mot.clust <- as.data.table(readxl::read_xlsx("db/motif/motif_annotations.xlsx", sheet= 1))
mot.id <- as.data.table(readxl::read_xlsx("db/motif/motif_annotations.xlsx", sheet= 2))
clust <- merge(mot.clust[, .(ID= as.character(Cluster_ID), cluster= Name)], 
               mot.id[, .(ID= as.character(Cluster_ID), motif= Motif)])

# Import row labels ----
lab <- readRDS(seq.info)[, line.idx:= .I]
lab[, label:= factor(label, label.lvls)]

# COmpute motif counts ----
# pwms <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/jeff_non-redundant_pwms_jaspar_all.rds")
# counts <- vl_motifCounts(sequences = lab$sequence,
#                          pwm_log_odds= pwms,
#                          p.cutoff = 1e-4,
#                          genome= "mm10",
#                          bg= "genome")
# saveRDS(counts,
#         "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/motifs/motif_counts_genome_bg.rds")
# pwms <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/jeff_non-redundant_pwms_jaspar_all.rds")
# counts <- vl_motifCounts(sequences = lab$sequence,
#                          pwm_log_odds= pwms,
#                          p.cutoff = 1e-4,
#                          genome= "mm10",
#                          bg= "subject")
# saveRDS(counts,
#         "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/motifs/motif_counts_subject_bg.rds")

# Import motif counts matrix ----
counts <- as.data.table(readRDS(mot.counts.file)/(lab$width/1000))

# Compute motif enrichment in heart ----
counts.list <- split(counts, lab$label)
enr <- vl_motifEnrich(counts= counts.list["vista_ts"],
                      control.counts = counts.list[["rdm"]])
# enr <- readRDS(enr.file)
clust <- merge(clust,
               enr[, .(motif= name, log2OR, padj)],
               by= "motif")

# Retrieve motifs contrib in heart ----
contrib <- as.data.table(readRDS(contrib.file), keep.rownames = "motif")
setnames(contrib, tissue, "contrib")
clust <- merge(clust,
               contrib[, .(motif, contrib)],
               by= "motif")
# Max per motif cluster
clust[, contrib:= max(contrib), cluster]

# Select best motif per cluster (highest enrichment vs random) ----
clust[, log10padj:= ifelse(log2OR>0, -log10(padj), log10(padj))]
setorderv(clust, "log10padj", -1)
sel <- clust[contrib>contrib.cutoffs[2] | (log2OR>1 & padj<padj.cutoffs[2]), {
  .SD[1]
}, cluster]

# Check intersection between vista enrich and ----
sel[, vista:= padj<padj.cutoffs[1] & log2OR>log2(1.5)]
sel[, high.contrib:= contrib > contrib.cutoffs[1]]
sel[, class:= fcase(vista & high.contrib, "High contrib. + Enr. vista (vs rdm genomic)",
                    vista, "Enr. vista (vs rdm genomic)",
                    high.contrib, "High contrib.")]
sel[, class:= paste0(class, " (n= ", .N, ")"), class]

# Subset relevant motif counts ----
counts <- counts[, sel$motif, with= FALSE]
colnames(counts) <- sel$cluster

# Measure enrichment vs random genomic sequences subgroup ----
counts.list <- split(counts, 
                     lab$label)
mot.enr <- list()
for(ctl.group in c("dinuc.adj", "rdm", "vista_inactive_all_tissues")) {
  # Enrich
  ctl <- which(names(counts.list)==ctl.group)
  enr <- vl_motifEnrich(counts = counts.list[-ctl],
                        control.counts = counts.list[[ctl]])
  enr[, cl:= factor(cl, intersect(label.lvls, cl))]
  # OR
  log2OR <- dcast(enr,
                  cl~motif,
                  value.var = "log2OR")
  log2OR <- as.matrix(log2OR, 1)
  # Padj
  padj <- dcast(enr,
                cl~motif,
                value.var = "padj")
  padj <- as.matrix(padj, 1)
  stars <- apply(padj, 2, function(x) {
    ifelse(x<5e-2, "*", "")
  })
  # Store
  mot.enr[[ctl.group]] <- list(log2OR= log2OR,
                               padj= padj,
                               stars= stars)
}

# Subset sequences (N= tissue spe. vista enhancers) ----
Nseq <- sum(lab[[label.column]]=="vista_ts")
lab[, line.sel:= seq(.N) <= Nseq, label]
counts <- counts[lab$line.sel, ]
lab <- lab[lab$line.sel, ]

# Subset 25 sequences per label for hclust ----
lab[, line.sel:= seq(.N)<=25, label] 
sub <- counts[lab$line.sel, ]
sub.lab <- lab[lab$line.sel, ]

# Plots ----
leg.cex <- .6
# Plot contrib vs enrich ----
pdf("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/pdf/motif_enriched_designed_seq.pdf",
    width = 7.75,
    height= 5.75)
vl_par(mfrow= c(2,2),
       xpd= F)
# Contrib vs enrichment
plot(clust$contrib,
     ifelse(clust$log2OR>0, -log10(clust$padj), log10(clust$padj)),
     xlab= "Contrib",
     ylab= "p.adjust (log10)",
     pch= 16,
     cex= .7)
abline(v= contrib.cutoffs,
       col= c("limegreen", "red"))
abline(h= -log10(padj.cutoffs),
       col= c("limegreen", "red"))
plot(sort(clust$contrib),
     ylab= "Contrib",
     pch= 16,
     cex= .7)
abline(h= contrib.cutoffs,
       col= c("limegreen", "red"))
plot(sort(ifelse(clust$log2OR>0, -log10(clust$padj), log10(clust$padj))),
     ylab= "Enrichment p.adj (-log10)",
     pch= 16,
     cex= .7)
abline(h= -log10(padj.cutoffs),
       col= c("limegreen", "red"))

# Overlaps
vl_par(mai= c(1.5, 2,.9,.9))
upsetPlot(list(high.contrib= sel[(high.contrib), motif],
               enriched.vista= sel[(vista), motif]))

# Heatmaps
layout(matrix(1:4, ncol= 1),
       heights = c(3.5,1,1,1/4.5))
vl_par(mai= c(.02, 1.75, .02, 2),
       omi= c(1,0,1,0),
       xpd= NA)
col <- vl_heatmap(counts,
                  cluster.rows = lab$label,
                  show.row.clusters = "left",
                  col= c("blue", "yellow", "red"),
                  breaks = seq(0, 6),
                  useRaster= T,
                  cluster.cols = T,
                  show.colnames = F,
                  legend.title = "Motif counts/kb",
                  legend.cex = leg.cex,
                  main= "heart",
                  col.annotations = sel$class,
                  gap.width = 1/100)$cols
idx <- col[(order), name]

# Enrichments
vl_heatmap(mot.enr$rdm$log2OR[, idx],
           show.numbers = mot.enr$rdm$stars[, idx],
           cluster.rows = F,
           show.colnames = F,
           legend.title = "OR vs. rdm genomic seq. (log2)",
           legend.cex = leg.cex)
vl_heatmap(mot.enr$dinuc.adj$log2OR[, idx],
           show.numbers = mot.enr$dinuc.adj$stars[, idx],
           cluster.rows = F,
           show.colnames = F,
           legend.title = "OR vs. rdm starting seq. (dinuc. adj.) (log2)",
           legend.cex = leg.cex)
# vl_heatmap(mot.enr$vista_inactive_all_tissues$log2OR[, idx],
#            show.numbers = mot.enr$vista_inactive_all_tissues$stars[, idx],
#            cluster.rows = F,
#            show.colnames = F,
#            legend.title = "OR vs. vista inact. all tissues (log2)",
#            legend.cex = leg.cex)

# Contribution
contribs <- matrix(sel$contrib, nrow= 1)
colnames(contribs) <- sel$cluster
rownames(contribs) <- "Contribution"
vl_heatmap(contribs[, idx, drop= F],
           legend.title = "Contrib",
           legend.cex = leg.cex,
           breaks= seq(-0.001, 0.001, length.out= 21))

# Clustering
vl_par(mfrow= c(1,1),
       mai= c(.9, 1.75, 1, 2),
       omi= c(0,0,0,0))
vl_heatmap(sub[, idx, with= F],
           cluster.rows = T,
           show.row.clusters = "left",
           col= c("blue", "yellow", "red"),
           breaks = seq(0, 6),
           useRaster= T,
           legend.title = "Motif counts/kb",
           legend.cex = leg.cex,
           main= "heart",
           row.annotations = sub.lab$label,
           row.annotations.col = label_colors,
           cutree.rows = 20,
           gap.width = 1/100)
vl_heatmap(sub[, idx, with= F][, sel[contrib>0.00025, cluster], with= F],
           cluster.rows = T,
           kmeans.k = 10,
           show.row.clusters = "left",
           col= c("blue", "yellow", "red"),
           breaks = seq(0, 6),
           useRaster= T,
           legend.title = "Motif counts/kb",
           legend.cex = leg.cex,
           main= "heart",
           row.annotations = sub.lab$label,
           row.annotations.col = label_colors)
dev.off()
