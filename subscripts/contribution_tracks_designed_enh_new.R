setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require("BSgenome.Mmusculus.UCSC.mm10")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Metadata ----
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

# Import contributions scores ----
dat <- meta[, {
  importContrib(
    contrib.file,
    fa = fa.file,
    FUN = function(x) colSums(x)
  )
}, tissue]
sel <- dat[(tissue=="heart" & seqnames=="311") | (tissue=="limb" & seqnames=="51") | (tissue=="midbrain" & seqnames=="5")]
setnames(sel, "tissue", "seqlvls")

# Call motif positions ----
TF <- readRDS("Rdata/motif_clusters_paper_3_tissues.rds")
TF <- TF[, .(TF= .(unique(unlist(TF)))), .(name, cluster)]
mot.file <- "db/motifs/motifs_validated_sequences_paper.rds"
if(!file.exists(mot.file)) {
  all_TFs <- readRDS("/groups/stark/vloubiere/motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")
  seq <- sel$seq
  names(seq) <- sel$seqlvls
  mot <- vl_motifPos(sequences = seq,
                     pwm_log_odds = all_TFs$pwms_log_odds[match(TF$name, all_TFs$meta$ID)],
                     bg = "genome",
                     genome = "mm10",
                     p.cutoff = 5e-4)
  saveRDS(mot, mot.file)
} else
  mot <- readRDS(mot.file)

# Simplify motifs
mot <- mot[!is.na(mot.count)]
mot <- mot[, ir[[1]], .(seqnames= seqlvls, motif)]
mot <- merge(mot, TF[, .(name, cluster)], by.x= "motif", by.y= "name")
mot <- mot[
  (seqnames=="heart" & cluster==1) | (seqnames=="limb" & cluster==2) |
    (seqnames=="midbrain" & cluster==3) | cluster=="4"
]

# Filter motifs overlapping seqlets ----
seqlets <- contribSeqlets(sel, zscore.cutoff = .25, log2OR.cutoff = 0, FDR.cutoff = .1)
setnames(seqlets, "seqlvls", "seqnames")
sel.mot <- intersectBed(mot, seqlets)

# Collapse per TF name ----
sel.mot[, TF:= tstrsplit(motif, "_", keep= 3)]
sel.mot[, TF:= toupper(TF)]
sel.mot <- sel.mot[, collapseBed(.SD), TF]

# Plot ----
sel[, idx:= .I]
# sel[, screen.start:= c(75, 350, 207)]
# sel[, screen.end:= screen.start+250]
sel[, screen.start:= 1]
sel[, screen.end:= screen.start+1000]

sel[, screen.start:= c(85, 550, 225)]
sel[, screen.end:= screen.start+175]

pdf("pdf/0_paper/activity_contribution_tracks_validated_enh_2.pdf", width = 10, height = 9)
vl_par(mfrow= c(3,1))
sel[, {
  
  # Plot contrib
  contribSeqLogo(contrib = .SD, xaxt= "n", start = screen.start, end= screen.end)
  title(main= seqlvls)
  
  # subset motifs
  .c <- sel.mot[.BY, on= "seqnames==seqlvls"]
  .c[, col:= adjustcolor(rainbow(.NGRP))[.GRP], TF]
  .c[, adjust:= seq(.N), collapseBed(.c, return.idx.only = T)]
  
  # Plot motifs
  top <- min(contrib.score[[1]])
  rect(
    xleft = .c$start-1,
    ybottom = top-strheight("M")*.c$adjust,
    xright = .c$end,
    ytop = top-strheight("M")*(.c$adjust-1),
    col= .c$col,
    border= NA,
    xpd= T
  )
  abline(v= seq(100, 1000, 100))
  # Plot TF name
  text(
    x = .c[, start+end]/2,
    y= top-strheight("M")*(.c$adjust-0.5), 
    .c$TF,
    xpd= T
  )
  
  # subset seqlets
  # .s <- seqlets[.BY, on= "seqnames==seqlvls"]
  # rect(
  #   xleft = .s$start-1,
  #   ybottom = top,
  #   xright = .s$end,
  #   ytop = top-strheight("M")*2,
  #   xpd= T
  # )
}, seqlvls]
dev.off()
