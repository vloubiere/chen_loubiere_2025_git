setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Selected sequences object were created in:
# file.edit("git_deepATAC/subscripts/create_clean_list_sequences.R") # Clean list of designed sequences (different EVO/LEDIDI designs...)
# file.edit("git_deepATAC/subscripts/set_cutoffs_predicted_activity_evgeny.R") # BLAST/specificity filters for LEDIDI_12_14...

# Import objects ----
dat <- list(
  heart= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_heart.rds"),
  limb= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_limb.rds"),
  midbrain= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_midbrain.rds")
)
dat <- rbindlist(dat, idcol = "tissue")
dat <- dat[label=="ledidi_12_14"] # This is the approach we used

# Fix labels ----
dat[, tissue:= factor(tissue, c("heart", "limb", "midbrain"))]
setorderv(dat, "tissue")
setnames(dat,
         c("act_heart", "act_limb", "act_midbrain"),
         c("Heart", "Limb", "Midbrain"))

# Plot ----
# pdf("pdf/0_paper/boxplot_compare_predicted_act_per_tissue.pdf", 3.1, 2.8)
pdf("pdf/_revision/boxplot_compare_predicted_act_per_tissue.pdf", 3.1, 2.8)
vl_par(mfrow= c(1,3),
       mai= c(.9, .4, .9, .1),
       lwd= .5)
dat[, {
  vl_boxplot(
    .SD,
    outline= T,
    ylab= "Pred. activity score (logits)",
    violin= T,
    main= paste0("Target tissue:\n", tissue),
    names= gsub("act_(.*)", "", names(.SD)),
    tilt.names= T,
    pts.cex= .1,
    whisklty= 0,
    viowex= .5,
    boxwex= .15
  )
  if(tissue=="heart") {
    abline(h= c(0,6), lty= 3)
    sel.enh <- unlist(.SD[id %in% c(311, 726, 834, 890, 845)])
  }
  if(tissue=="limb") {
    abline(h= c(0,5), lty= 3)
    sel.enh <- unlist(.SD[id %in% c(1104,  1121,  51,  103,  112)])
  }
  if(tissue=="midbrain") {
    abline(h= c(0,7), lty= 3)
    sel.enh <- unlist(.SD[id %in% c(1008, 734, 543, 169, 5)])
  }
  points(rep(1:3, each= 5)+.4, sel.enh, pch= 16, cex= .5, col= adjustcolor("limegreen", .7))
  print("")
}, tissue, .SDcols= c("Heart", "Limb", "Midbrain")]
dev.off()