setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# VISTA enhancers size distribution ----
# Import vista tiles
vista <- readxl::read_excel("/groups/stark/shenzhi.chen/db/VISTA_enhancer_dataset/VISTA2024_AllTissuesReferenceAlleles.xlsx")
vista <- as.data.table(vista)
vista <- vista[, .(class= "vistaTile",
                   peakID= Vista.ID,
                   genome= fcase(grepl("^hs", Vista.ID), "hg38",
                                 grepl("^mm", Vista.ID), "mm10"),
                   coor_hg38= Coordinates_hg38,
                   coor_mm10= Coordinates.mm10,
                   heart= Heart,
                   limb= Limb,
                   forebrain= Forebrain,
                   midbrain= Midbrain,
                   hindbrain= Hindbrain,
                   neuralTube= NeuralTube)]
vista <- vista[!is.na(genome)]
# Remove mutated enhancers
mut <- fread("db/public/VISTA_experiments_20250609.tsv")[grepl("mutagenesis", description), vista_id]
vista <- vista[!(peakID %in% mut)]
# Extract genomic coordinates
vista[genome=="mm10", c("seqnames", "start", "end"):= importBed(coor_mm10)[, 1:3]]
vista[genome=="hg38", c("seqnames", "start", "end"):= importBed(coor_hg38)[, 1:3]]
# Width
vista[, width:= end-start+1]

# Compute minimum distance between ledidi edits ----
heart <- readRDS("Rdata/final_designed_enhancer_sequences_heart.rds")
heart <- heart[id %in% c(311, 726, 834, 890, 845) & grepl("ledidi_12_14", label)]
limb <- readRDS("Rdata/final_designed_enhancer_sequences_limb.rds")
limb <- limb[id %in% c(1104,  1121,  51,  103,  112) & grepl("ledidi_12_14", label)]
CNS <- readRDS("Rdata/final_designed_enhancer_sequences_midbrain.rds")
CNS <- CNS[id %in% c(1008, 734, 543, 169, 5) & grepl("ledidi_12_14", label)]
dat <- rbindlist(list(heart= heart, limb= limb, CNS= CNS), idcol = "tissue")
dat[, label:= switch(
  label,
  "ledidi_12_14"= "validated.enhancer",
  "ledidi_12_14_ini"= "initialization.seq"
), label]
dat <- dat[, .(tissue, label, init.id= id, sequence)]
dat <- dcast(dat, tissue+init.id~label, value.var = "sequence")
dat[, c("min", "max"):= {
  V1 <- strsplit(initialization.seq, "")[[1]]
  V2 <- strsplit(validated.enhancer, "")[[1]]
  .(min(which(V1!=V2)), max(which(V1!=V2)))
}, .(initialization.seq, validated.enhancer)]
dat[, width:= max-min+1]
edit <- dat[, .(width= min(width)), tissue]
edit[, at:= switch(tissue, "heart"= 2, "CNS"= 4, "limb"= 3), tissue]

# Compute minimum distance between important seqlets ----
# Metadata
meta <- data.table(tissue= c("heart", "limb", "midbrain"))
meta[, fa.file:= paste0(
  "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/Leididi_design/",
  tissue,
  "/selected_sequences/all_ledidi_design_sequence.fasta"
)]
meta[, contrib.file:= paste0(
  "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/Leididi_design/",
  tissue,
  "/selected_sequences/act_contri/Model_final_selected_ledidi_designed_sequence.fasta_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5"
)]
# Import contrib
contrib <- meta[, importContrib(h5 = contrib.file, fa= fa.file, FUN = function(x) colSums(x)), tissue]
contrib[, id:= as.integer(seqnames)]
contrib <- contrib[
  (tissue=="heart" & id %in% c(311, 726, 834, 890, 845)) |
    (tissue=="limb" & id %in% c(1104,  1121,  51,  103,  112)) |
    (tissue=="midbrain" & id %in% c(1008, 734, 543, 169, 5))
]
contrib[, c("min", "max"):= {
  .c <- scale(contrib.score[[1]])
  .(min(which(.c>1)), max(which(.c>1)))
}, .(tissue, seqnames)]
contrib[, width:= max-min+1]
contrib <- contrib[, .(width= min(width)), tissue]
contrib[, at:= switch(tissue, "heart"= 2, "midbrain"= 4, "limb"= 3), tissue]

# Plot ----
pdf("pdf/_revision/vista_size_distribution.pdf", 4.25, 3.5)
vl_par(mai= c(.9, .9, .9, 1.8+.5))
vl_boxplot(
  list(
    All= vista[, width],
    Heart= vista[heart==1, width],
    Limb= vista[limb==1, width],
    Midbrain= vista[midbrain==1, width]
  ),
  main= "VISTA sequences",
  violin = T,
  tilt.names = T,
  ylab= "Genomic size (bp)",
  boxwex= .2,
  viowex= .8,
  viocol= "lightgrey",
  col= "white"
)
abline(h= 1001, lty= "11")
text(par("usr")[2], 1001, 1001, xpd= T, pos= 4, cex= .5, offset= 0.2)
abline(h= min(edit$width), lty= "11", col= "red")
text(par("usr")[2], min(edit$width), min(edit$width), xpd= T, pos= 4, cex= .5, offset= 0.2, col= "red")
edit[, points(at, width, pch= 16, col= "red", cex= .5)]
contrib[, points(at, width, pch= 16, col= "limegreen", cex= .5)]
vl_legend(
  legend= c("Min. dist. between edits", "Min. dist. between important seqlets (zscore>1)"),
  pch= 16,
  col= c("red", "limegreen"),
  title= "Validated synthetic enhancers"
)
dev.off()