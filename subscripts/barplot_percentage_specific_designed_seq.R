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
dat[, clean.tissue:= switch(as.character(tissue), "heart"= "Heart", "limb"= "Limb", "midbrain"= "Midbrain (CNS)"), tissue]
perc <- dat[(active), .(perc= sum(specific)/.N*100), clean.tissue]

# Plot ----
pdf("pdf/0_paper/barplot_percentage_specific_designed_seq.pdf", 1.8, 2.8)
vl_par(mai= c(.9, .9, .9, .4),
       lwd= .5)
bar <- vl_barplot(perc$perc,
                  ylab= "% pred. tissue-specific\n(<0 other tissues)",
                  names.arg= perc$clean.tissue)
text(bar$x, bar$y, round(perc$perc, 1), xpd= T, cex= 5/12, pos= 3, offset= .25)
dev.off()