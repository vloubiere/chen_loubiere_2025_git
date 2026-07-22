setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import validated synthtetic enhancer sequences ----
heart <- readRDS("Rdata/final_designed_enhancer_sequences_heart.rds")
heart <- heart[id %in% c(311, 726, 834, 890, 845) & label=="ledidi_12_14"]
limb <- readRDS("Rdata/final_designed_enhancer_sequences_limb.rds")
limb <- limb[id %in% c(1104,  1121,  51,  103,  112) & label=="ledidi_12_14"]
CNS <- readRDS("Rdata/final_designed_enhancer_sequences_midbrain.rds")
CNS <- CNS[id %in% c(1008, 734, 543, 169, 5) & label=="ledidi_12_14"]
dat <- rbindlist(list(heart= heart, limb= limb, CNS= CNS), idcol = "tissue")
dat <- dat[, .(tissue, id, sequence)]

# Count motifs ----
mot.db <- readRDS("../../motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")
mot <- vl_motifCounts(
  sequences = dat$sequence,
  pwm_log_odds = mot.db$pwms_log_odds,
  # bg= "subject",
  bg= "genome",
  genome = "mm10",
  p.cutoff = 1e-4
)

# Make table ----
final <- cbind(dat[, !"sequence"], as.data.table(mot))
final <- melt(final, id.vars = c("tissue", "id"))

# Save ----
saveRDS(final, "db/motifs/revision_motif_counts_validated_synthetic_enhancer_seq.rds")