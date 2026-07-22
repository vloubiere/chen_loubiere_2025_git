setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
folder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/"
clust <- readRDS(paste0(folder, "Rdata/final_enhancer_selection/all_motifs_clusters_metatable.rds"))
# Remove Heart contrib scores
clust$contrib <- clust$access <- NULL

# Rename relevant motif clusters ----
# Accessibility
clust[grepl("CREB", cluster, ignore.case = T), new.name:= paste0("CREB.accessibility")]
clust[grepl("AP1", cluster, ignore.case = T), new.name:= paste0("AP1.accessibility")]
clust[grepl("ZIC", cluster, ignore.case = T), new.name:= paste0("ZIC.accessibility")]
clust[grepl("KLF/SP/", cluster, ignore.case = T), new.name:= paste0("KLF/SP.accessibility")]
# Heart
clust[grepl("NKX2.5", motif, ignore.case = T), new.name:= "NKX2-5.heart"]
clust[grepl("HOXB2|HOXB3|HOXB5", motif, ignore.case = T), new.name:= "HOXB2-5.heart"]
clust[grepl("IRX", motif, ignore.case = T), new.name:= "IRX.heart"]
clust[grepl("TBX4", motif, ignore.case = T), new.name:= "TBX4.heart_limb"]
clust[cluster %in% c("TFAP2/1", "TFAP2/2"), new.name:= "TFAP2.heart"]
clust[cluster %in% c("GATA", "MEF2", "TEAD", "EGR", "NRF1", "HAND1", "Ebox/CACGTG/1", "Ebox/CACGTG/2"), new.name:= paste0(cluster, ".heart")]
clust[cluster %in% c("PAX-halfsite", "PRDM16"), new.name:= paste0(cluster, ".heart_unexpected")]
# Limb
clust[grepl("HOXD9|HOXD10|HOXD11|HOXD13", motif, ignore.case = T), new.name:= "HOXD9-13.limb"]
clust[grepl("GLI3", motif, ignore.case = T), new.name:= "GLI3.limb"]
clust[grepl("LMX1B", motif, ignore.case = T), new.name:= "LMX1B.limb"]
clust[grepl("SOX9", motif, ignore.case = T), new.name:= "SOX9.limb"]
clust[cluster %in% c("RUNX/1","RUNX/2"), new.name:= "RUNX.limb"]
clust[cluster %in% c("Ebox/CAGCTG"), new.name:= "TWIST1.limb"]
clust[grepl("olig3", motif, ignore.case = T), new.name:= "olig3.limb"]
clust[grepl("HTF4", motif, ignore.case = T), new.name:= "HTF4.limb"]
# Midbrain
clust[grepl("OTX2", motif, ignore.case = T), new.name:= "OTX2.midbrain"]
clust[grepl("SOX9", motif, ignore.case = T), new.name:= "SOX9.midbrain"]
clust[grepl("OTX1", motif, ignore.case = T), new.name:= "OTX1.midbrain"]
clust[grepl("RFX", cluster, ignore.case = T), new.name:= "RFX.midbrain"]
clust[grepl("POU/1", cluster, ignore.case = T), new.name:= "POU.midbrain"]
clust[grepl("nkx21", motif, ignore.case = T), new.name:= "NKX21.midbrain"]
clust[grepl("NEUROD1", motif, ignore.case = T), new.name:= "NEUROD1.midbrain"]
clust[grepl("OLIG2", motif, ignore.case = T), new.name:= "OLIG2.midbrain"]
clust[grepl("Pax6", motif, ignore.case = T), new.name:= "PAX6.midbrain"]
clust[grepl("FOXG1", motif, ignore.case = T), new.name:= "FOXG1.midbrain"]
clust[grepl("Ngn2", motif, ignore.case = T), new.name:= "NGN2.midbrain"]
clust[grepl("TBR1_TBX_2", motif, ignore.case = T), new.name:= "TBR1.midbrain"]
clust[grepl("FEZF1", motif, ignore.case = T), new.name:= "FEZF.midbrain"]
clust[is.na(new.name), new.name:= paste0(cluster, ".other")]
clust <- clust[!is.na(new.name)]

# Import sequences ----
seq <- list(heart= "/Rdata/final_enhancer_selection/final_designed_enhancer_sequences_heart.rds",
            limb= "/Rdata/final_enhancer_selection/final_designed_enhancer_sequences_limb.rds",
            midbrain= "/Rdata/final_enhancer_selection/final_designed_enhancer_sequences_midbrain.rds")
seq <- lapply(seq, function(x) readRDS(paste0(folder, x)))
seq <- lapply(seq, function(x) {
  x[label %in% c("ledidi_12_14_ini", "ledidi_12_14"), .(label, selected, id, sequence)]
})
seq <- rbindlist(seq, idcol = "tissue")

# Compute motif counts ----
mot <- readRDS(paste0(folder, "db/motif/jeff_non-redundant_pwms_jaspar_all.rds"))
mot <- mot[sapply(mot, name) %in% clust$motif]
counts.file <- "db/motifs/motif_counts_designed_sequences_0.0001.rds"
if(!file.exists(counts.file)) {
  counts <- vl_motifCounts(sequences = seq$sequence,
                           pwm_log_odds = mot,
                           genome = "mm10",
                           bg = "subject",
                           p.cutoff = 1e-4)
  saveRDS(counts, counts.file)
}
counts <- readRDS(counts.file)

# Subtract initiation sequence motif counts to designed sequences  ----
diff <- counts[seq$label=="ledidi_12_14",]-counts[seq$label=="ledidi_12_14_ini",]
diff[, tissue:= seq[label=="ledidi_12_14", tissue]]
diff[, selected:= seq[label=="ledidi_12_14", selected]]
diff <- melt(diff,
             id.vars = c("tissue", "selected"),
             variable.name = "motif")

# Add cluster info ----
dat <- merge(diff,
              clust,
              by= "motif")

# Compute stats  ----
dat <- dat[, .(mean= mean(value), perc= sum(value>0)/.N*100), .(new.name, motif, tissue, selected)]
dat[, sel.mean:= motif %in% dat[mean>.5, motif[which.max(mean)], new.name]$V1]
dat[, sel.perc:= motif %in% dat[perc>40, motif[which.max(mean)], new.name]$V1]

# Mean matrix ----
mat.mean <- dcast(dat[(sel.mean)],
                  new.name~tissue+selected,
                  value.var = "mean")
mat.mean <- as.matrix(mat.mean, 1)
annot.mean <- unlist(tstrsplit(rownames(mat.mean), "[.]", keep= 2))
rownames(mat.mean) <- unlist(tstrsplit(rownames(mat.mean), "[.]", keep= 1))

# Perc matrix ----
mat.perc <- dcast(dat[(sel.perc)],
                  new.name~tissue+selected,
                  value.var = "perc")
mat.perc <- as.matrix(mat.perc, 1)
annot.perc <- unlist(tstrsplit(rownames(mat.perc), "[.]", keep= 2))
rownames(mat.perc) <- unlist(tstrsplit(rownames(mat.perc), "[.]", keep= 1))

# Plot heatmaps ----
pdf("pdf/frequence_all_known_motifs.pdf", width = 7.75, height = 8.5)
vl_par(mai= c(.8, .9, .9, 1.9),
       mfrow= c(1,2))
br <- seq(-1, 4, .1)
col <- circlize::colorRamp2(c(min(br), 0, max(br)), colors = c("cornflowerblue", "white", "tomato"))(br)
vl_heatmap(mat.mean,
           breaks = br,
           col = col,
           show.numbers = round(mat.mean, 1),
           row.annotations = annot.mean,
           legend.title = "Mean counts",
           numbers.cex = .5)
br <- seq(0, 80, 10)
col <- circlize::colorRamp2(range(br), colors = c("white", "red"))(br)
vl_heatmap(mat.perc,
           breaks = br,
           col = col,
           show.numbers = round(mat.perc),
           row.annotations = annot.perc,
           legend.title = "% seq >=1 motif",
           numbers.cex = .5)
dev.off()
