setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
folder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/"
clust <- readRDS(paste0(folder, "Rdata/final_enhancer_selection/all_motifs_clusters_metatable.rds"))

# Rename HD cluster... ----
clust[grepl("NKX2.5", motif, ignore.case = T), cluster]
clust[cluster %in% c("HD/4", "HD/22"), cluster:= "NKX2.5"]
clust[grepl("HOXB2|HOXB3|HOXB5", motif, ignore.case = T)]
clust[cluster %in% c("HD/2"), cluster:= "HOXB2-5"]
clust[grepl("HOXD9|HOXD10|HOXD11|HOXD13", motif, ignore.case = T)]
clust[cluster %in% c("HD/17", "HD/18"), cluster:= "HOXD9-13"]
clust[grepl("HOXA9|HOXA10|HOXA11|HOXA13", motif, ignore.case = T)]
clust[cluster %in% c("HD/14"), cluster:= "HOXA9"]
clust[grepl("IRX", motif, ignore.case = T)]
clust[cluster %in% c("HD/23"), cluster:= "IRX"]

# Heart motifs ----
heart <- list(
  `Heart-specific`= c("GATA", "MEF2", "TEAD", "ZIC", "EGR", "NRF1"),
  General= c(unique(grep("Ebox|CREB|AP1", clust$cluster, value = T)), "KLF/SP/1", "KLF/SP/2", "TFAP2/1", "TFAP2/2"),
  Missing= c(unique(grep("TBX/|FOX/", clust$cluster, value = T)), "HAND1", "NKX2.5", "HOXB2-5", "IRX"),
  Unexpected= c("PAX-halfsite", "PRDM16"),
  Limb= c("HOXD9-13", "HOXA9")
)
heart <- rbindlist(lapply(heart, as.data.table), idcol = "Class")
setnames(heart, "V1", "cluster")
heart <- merge(heart,
               clust,
               all.x= T,
               by= "cluster")
if(any(is.na(heart$motif)))
  stop("Motifs did not match/missing")

# Import enrichments ----
enr <- readRDS("db/motifs/enrichment_selected_sequences.rds")
enr <- enr[motif %in% heart$motif & cl %in% c("vista_ts", "ledidi_12_14")]
enr[, cl:= droplevels(cl)]
enr[heart, name:= i.cluster, on= "motif"]
enr <- enr[motif %in% enr[cl=="vista_ts", motif[which.min(padj)], name]$V1]

# Enrichment vs random sequences ----
counts <- readRDS("db/motifs/motif_counts_selected_sequences_0.0001.rds")
enr2.file <- "db/motifs/enrichment_ledidi_vs_ini.rds"
if(!file.exists(enr2.file)) {
  enr2 <- vl_motifEnrich(counts[label=="ledidi_12_14", -c("label", "id", "width", "selected")],
                         counts[label=="ledidi_12_14_ini", -c("label", "id", "width", "selected")])
  saveRDS(enr2, enr2.file)
}
enr2 <- readRDS(enr2.file)
enr2 <- enr2[motif %in% heart$motif]
enr2[heart, name:= i.cluster, on= "motif"]
enr2 <- enr2[, {
  if(any(padj<0.05 & log2OR>0))
    .SD[log2OR>0][which.min(padj)] else
      .SD[which.min(padj)]
}, name]
enr2[, padj:= ifelse(padj==0, min(padj[padj>0]), padj)]

# Plot
pdf("pdf/enrichment_known_heart_motifs.pdf", width = 4, height = 6)
vl_par(mai= c(.9, 1.5, .9, 1.5))
plot(enr,
     padj.cutoff = .05,
     top.enrich = Inf,
     cex = .4,
     plot.empty.clusters = T)
plot(enr2,
     log2OR.abs.cutoff = 0,
     padj.cutoff = 1,
     top.enrich = Inf,
     breaks = seq(0, 10),
     col= c("pink", "red"))
dev.off()
