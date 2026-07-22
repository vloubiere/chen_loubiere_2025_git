setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import contributions
contrib <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/motifs_enrichment_analysis/ledidi_act_acc/motifs_collapsed_heart_mtx.rds")$contri_score

# Import motifs matrix ----
mot <- readRDS("db/motifs/individual_motif_counts_selected_sequences.rds")
setnames(mot, c("contrib", "access"), c("enh.contrib", "acc.contrib"))

# Retrieve labels ----
labels <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/motifs_enrichment_analysis/ledidi_act_acc/motifs_collapsed_heart_mtx.rds")$counts_mtx
mot[labels, label:= i.label, on= "seqnames"]

# Remove useless labels ----
mot <- mot[label!="vista_other_tissues"]
mot[, label:= factor(label,
                     c("dinuc.adj", "rdm",  "vista_inactive_all_tissues", "vista_ts", "ledidi_design(12,10)", "ledidi_selected"))]

# Select motifs with high score ----
enh.contrib.cutoffs <- 0.0004
acc.contrib.cutoffs <- 0.0015
mot <- mot[enh.contrib>enh.contrib.cutoffs | acc.contrib>acc.contrib.cutoffs]

# Select motifs using enrichment ----
sel.file <- "db/motifs/select_motifs_with_strongest_enrichment_vista_vs_rdm.rds"
if(!file.exists(sel.file)) {
  sel <- dcast(mot[label %in% c("vista_ts", "rdm")],
               label+seqnames~motif,
               value.var = "score",
               fun.aggregate = length)
  sel$seqnames <- NULL
  sel <- split(sel[, !"label"], sel$label, drop = TRUE)
  enr <- vl_motifEnrich(sel$vista_ts,
                        sel$rdm)
  enr[unique(mot[, .(cluster, motif)]), cluster:= i.cluster, on= "motif"]
  saveRDS(enr, sel.file)
}
sel.mot <- readRDS(sel.file)
sel.mot <- sel.mot[, .SD[which.min(padj)], cluster]
mot <- mot[motif %in% sel.mot$motif]
mot <- dcast(mot,
             seqnames+label~cluster,
             value.var = "score",
             fun.aggregate = length)
mot$seqnames <- NULL
setorderv(mot, "label")

# Compute motif enrichment ----
enr.file <- "db/motifs/enrichment_selected_sequences.rds"
if(!file.exists(enr.file)) {
  # Remove useless clusters
  sub <- mot[label!="ledidi_selected"]
  # Split
  mot.list <- split(sub[, !"label"],
                    sub[, label],
                    drop = TRUE)
  # Compute enrichment
  enr <- vl_motifEnrich(counts = mot.list[names(mot.list) != "rdm"],
                        control.counts = mot.list[["rdm"]])
  saveRDS(enr, enr.file)
}
enr <- readRDS(enr.file)
# Order based on enrichment
enr[, name:= factor(name, enr[cl=="vista_ts"][order(log2OR), name])]
enr[, cl:= factor(cl, 
                  c("dinuc.adj", "vista_inactive_all_tissues", "vista_ts", "ledidi_design(12,10)"))]

# Heatmap enrichments ----
OR <- dcast(enr,
            cl~name,
            value.var = "log2OR")
OR <- as.matrix(OR, 1)
pval <- dcast(enr,
              cl~name,
              value.var = "padj")
pval <- as.matrix(pval, 1)
stars <- apply(pval, 2, function(x) ifelse(x<0.005, "*", ""))

# Add contributions ----
mean.contrib <- contrib[cluster %in% enr$name, .(enh.contrib= mean(contrib), acc.contrib= mean(access)), cluster]
mean.contrib <- melt(mean.contrib, id.vars = "cluster")
mean.contrib[, cluster:= factor(cluster, levels(enr$name))]
contrib.mat <-  dcast(mean.contrib, variable~cluster, value.var = "value")
contrib.mat <- as.matrix(contrib.mat, 1)

# Add percentage of enhancer containing motif ----
perc <- mot[, lapply(.SD, function(x) sum(x>0)/.N*100), label]
perc <- as.matrix(perc, 1)
perc <- perc[, levels(enr$name)]
perc <- perc[c("dinuc.adj", "rdm", "vista_inactive_all_tissues", "vista_ts", "ledidi_design(12,10)"),]

# Plot ----
pdf("pdf/final_heatmap_selected_enhancers.pdf", width = 10, height = 10)
layout(mat= matrix(1:6, ncol= 1),
       heights= c(3.6,1.4,1.4,4.4,15,7))
vl_par(mai= c(0.1, 1.8, 0.1, 2),
       omi= c(.7, 0, .2, 0))
vl_heatmap(OR,
           cluster.rows = FALSE,
           legend.title = "log2OR vs.rdm genomic seq.",
           show.numbers = stars,
           show.colnames = FALSE)
vl_heatmap(contrib.mat[1,,drop= FALSE],
           cluster.rows = FALSE,
           legend.title = "Enh.contrib",
           show.colnames = FALSE,
           show.legend = FALSE)
vl_heatmap(contrib.mat[2,,drop= FALSE],
           cluster.rows = FALSE,
           legend.title = "Acc.contrib",
           show.colnames = FALSE,
           col= c("white", "red"),
           show.legend = FALSE)
vl_heatmap(perc,
           cluster.rows = FALSE,
           legend.title = "% seq. >=1 mot.",
           col= c("white", "red"),
           show.colnames = FALSE,
           show.numbers = round(perc))
# All enhancers
vl_heatmap(as.matrix(mot[label!="ledidi_selected", colnames(OR), with= FALSE], 1),
           cluster.rows = mot[label!="ledidi_selected", label],
           show.row.clusters = "left",
           breaks = seq(0, 10),
           col= colorRampPalette(c("blue", "yellow", "red"))(11),
           legend.title = "Motif counts",
           show.colnames = FALSE)
# Selected enhancers
vl_heatmap(as.matrix(mot[label=="ledidi_selected", colnames(OR), with= FALSE], 1),
           show.row.clusters = "left",
           breaks = seq(0, 10),
           col= colorRampPalette(c("blue", "yellow", "red"))(11),
           legend.title = "Motif counts")
dev.off()
