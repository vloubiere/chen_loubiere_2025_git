setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Initiate plot
pdf("pdf/selected_sequences_activity_evgeny.pdf", width = 30, height = 38)
vl_par(omi= c(.5,0,.5,1.5),
       mgp= c(.5, .25, 0))
layout(matrix(seq(12*3), ncol= 1),
       heights = c(1,1,.6))

# For each tissue ----
for(tiss in c("heart", "limb", "midbrain")) {
# for(tiss in "heart") {
  # Input files ----
  # fasta
  fa.file <- paste0(
    "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/Leididi_design/",
    tiss,
    "/selected_sequences/all_ledidi_design_sequence.fasta"
  )
  # Contribution
  enh.contrib.file <- paste0(
    "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/Leididi_design/",
    tiss,
    "/selected_sequences/acc_contri/Model_final_selected_ledidi_designed_sequence.fasta_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5"
  )
  acc.contrib.file <- paste0(
    "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/Leididi_design/",
    tiss,
    "/selected_sequences/act_contri/Model_final_selected_ledidi_designed_sequence.fasta_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5"
  )
  
  # Import enhancer contrib ----
  enh.contrib <- importContrib(h5 = enh.contrib.file,
                               # bed = tmp,
                               fa = fa.file,
                               FUN = function(x) colSums(x))
  
  # Import accessibility contrib ----
  acc.contrib <- importContrib(h5 = acc.contrib.file,
                               # bed = tmp,
                               fa = fa.file,
                               FUN = function(x) colSums(x))
  
  # Sequence info ----
  seq.info <- readRDS(paste0("Rdata/final_designed_enhancer_sequences_", tiss, ".rds"))
  designed.seq.info <- seq.info[label=="ledidi_12_14" & !is.na(selected_validation_name)][order(selected_validation_name)]
  ini.seq.info <- seq.info[label=="ledidi_12_14_ini"][designed.seq.info$id, on= "id"]
  
  # Subset to selected sequences ----
  enh.contrib <- enh.contrib[designed.seq.info$id, on= "seqnames"]
  enh.contrib$name <- designed.seq.info$selected_validation_name
  acc.contrib <- acc.contrib[designed.seq.info$id, on= "seqnames"]
  acc.contrib$name <- designed.seq.info$selected_validation_name
  print(paste(tiss, nrow(enh.contrib)))
  
  # Call positions selected motifs ----
  selected.motifs <- readRDS(paste0("Rdata/final_designed_enhancer_motifs_", tiss, ".rds"))
  mot <- readRDS("Rdata/annotated_PWMs.rds")
  mot <- mot[rev(unique(selected.motifs$motif)), pwm, on= "motif"]
  seq <- c(designed.seq.info$sequence, ini.seq.info$sequence)
  names(seq) <- c(designed.seq.info$id, paste0(ini.seq.info$id, "_ini"))
  mot.pos.file <- paste0("db/motifs/motif_pos_", tiss, "_sequences_0.0001_subject_bg.rds")
  if(!file.exists(mot.pos.file)) {
  # if(T) {
    mot.pos <- vl_motifPos(sequences = seq,
                           pwm_log_odds = mot,
                           bg = "subject",
                           p.cutoff = 1e-4,
                           scratch = "/scratch/stark/vloubiere/motifs/")
    saveRDS(mot.pos,
            mot.pos.file)
  }
  mot.pos <- readRDS(mot.pos.file)
  mot.pos <- vlite::motifPosToBed(mot.pos)
  
  # Split designed and initiation sequences ----
  ini.seq.mot <- mot.pos[grepl("_ini", seqlvls)]
  ini.seq.mot[, seqnames:= gsub("_ini$", "", seqlvls)]
  sel.seq.mot <- mot.pos[!grepl("_ini", seqlvls)]
  sel.seq.mot[, seqnames:= seqlvls]
  
  # Remove motifs that were present in initiation sequences ----
  # sel.seq.mot <- intersectBed(sel.seq.mot, ini.seq.mot, invert = T)
  
  # # Call seqlets ----
  # enh.seqLets <- vlite::contribSeqlets(enh.contrib, zscore.cutoff = 1, log2OR.cutoff = 0, FDR.cutoff = .1)
  # setnames(enh.seqLets, "seqlvls", "seqnames")
  # acc.seqLets <- vlite::contribSeqlets(acc.contrib, zscore.cutoff = 1, log2OR.cutoff = 0, FDR.cutoff = .1)
  # setnames(acc.seqLets, "seqlvls", "seqnames")
  
  # Remove motifs that do not overlap any significant seqlet ----
  # sel.seq.mot <- intersectBed(sel.seq.mot, rbind(enh.seqLets, acc.seqLets))
  
  # Only keep motifs that were there in the list of cluster characteristic motifs ----
  sel.seq.mot <- sel.seq.mot[motif %in% selected.motifs$motif]
  sel.seq.mot[selected.motifs, name:= i.cluster, on= "motif"]
  sel.seq.mot[, name:= factor(name, selected.motifs$cluster)]
  setorderv(sel.seq.mot, "name")
  
  # Add colors ----
  extra.colors <- c("black", "grey30", "grey70", "grey90",
                    "brown", "khaki4", "navajowhite4", "navajowhite2",
                    "darkslategray", "lightseagreen", "aquamarine3", "aquamarine",
                    "plum3", "plum2")
  Cc <- rainbow(length(levels(sel.seq.mot$name))-length(extra.colors)+1)
  Cc <- Cc[-length(Cc)]
  Cc <- c(Cc, extra.colors)
  sel.seq.mot[, col:= Cc[.GRP], name]
  
  # Collapse per simplified name ----
  sel.seq.mot <- sel.seq.mot[, {
    collapseBed(.SD)
  }, .(name, col), .SDcols= c("seqnames", "start", "end")]
  
  # Adjust y position for overlapping motifs ----
  sel.seq.mot$idx <- collapseBed(sel.seq.mot, return.idx.only = T)
  sel.seq.mot[, adjust:= seq(.N)-1, idx]
  
  # Plot ----
  for(i in seq(nrow(enh.contrib))) {
  # for(i in 1) {
    # Plot accessibility contrib
    par(mai= c(.01, 2.5, .01, .4))
    vlite::contribSeqLogo(contrib = acc.contrib,
                          row.idx = i,
                          xaxt= "n",
                          ylab= NA,
                          xlab= paste("Accessibility contrib.", acc.contrib[i, name]))
    # Labels and predicted act info
    mtext("Accessibility contrib.", side = 2, las = 1, line = 2)
    pred.act <- designed.seq.info[i, paste0("Predicted act: heart=", round(act_heart, 1), ", limb=", round(act_limb, 1), ", midbrain= ", round(act_midbrain, 1))]
    text(par("usr")[1], par("usr")[4], labels = paste0(enh.contrib[i, name], "; ", pred.act), pos= 4, cex= 1.5, xpd= NA)
    vl_legend("topright", legend= enh.contrib[i, name], cex= 1)
    # Title
    if(i==1)
      mtext(tiss, outer = T, cex= 1.5)
    # Plot enhancer contrib
    vlite::contribSeqLogo(contrib = enh.contrib,
                          row.idx = i,
                          xaxt= ifelse(i==nrow(enh.contrib), "s", "n"),
                          ylab= NA,
                          xlab= paste("Enhancer contrib.", enh.contrib[i, name]))
    mtext("Enhancer contrib.", side = 2, las = 1, line = 2)
    # Labels
    vl_legend("topright", legend= enh.contrib[i, name], cex= 1)
    # Plot Motifs
    height <- strheight("M", cex= 1.2)
    sel.seq.mot[enh.contrib$seqnames[i], {
      rect(xleft = start-1,
           ybottom = par("usr")[3]-height*(adjust+1),
           xright= end,
           ytop = par("usr")[3]-height*adjust,
           border= NA,
           col= adjustcolor(col, .6),
           xpd= NA)
      text(mean(c(start, end)),
           par("usr")[3]-height*(adjust+.5),
           name,
           cex= 1,
           xpd= NA)
    }, .(start, end, col, adjust, name), on= "seqnames"]
    # Legend
    if(i==1) {
      leg <- unique(sel.seq.mot[, .(col, name)])
      vl_legend(legend= leg$name,
                fill= adjustcolor(leg$col, .6),
                xpd= NA,
                cex= 1.5)
    }
    # Split line ----
    par(mai= c(0,0,0,0))
    plot.new()
    segments(0, 0.25, 1, .25)
  }
}
dev.off()