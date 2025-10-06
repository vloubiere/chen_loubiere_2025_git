setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# We found the following examples on IGV ----
coor <- c(
  "chr18:3,637,665-3,694,606", # Heart
  "chr18:62,159,828-62,186,459", # Heart + globally open
  "chr18:3,825,318-3,854,833", # Midbrain
  "chr18:60,915,859-60,935,423", # Midbrain
  "chr18:84,879,206-84,891,723", # Heart
  "chr18:69,599,055-69,635,656", # Limb
  "chr18:67,440,557-67,458,528", # globally open
  "chr18:65,797,620-65,803,440" # globally open
)
bed <- importBed(coor[c(7, 3, 1, 6)])

# Import PCC ----
PCC <- readRDS("db/PCC/obs_vs_pred_per_tissue_chr18_peaks_union.rds")
PCC[, PCC:= round(PCC, 2)]
setkeyv(PCC, "tissue")

# Plot ----
pdf("pdf/0_paper/screenshot_example_accessibility_models.pdf", width = 6, height = 3.5)
vl_par(mai= c(.9, 2, .9, .4))
vlite::bwScreenshot(
  bed = bed,
  tracks= c("db/bw/observed/midbrain_treat_pileup.bigwig",
            "db/bw/predicted/midbrain_predicted_accessibility_chr18.bw",
            "db/bw/observed/heart_treat_pileup.bigwig",
            "db/bw/predicted/heart_predicted_accessibility_chr18.bw",
            "db/bw/observed/limb_treat_pileup.bigwig",
            "db/bw/predicted/limb_predicted_accessibility_chr18.bw"),
  track.names = c(paste("PCC=", PCC["Midbrain", PCC], "| Midbrain obs."),
                  "Midbrain pred.",
                  paste("PCC=", PCC["Heart", PCC], "| Heart obs."),
                  "Heart pred.",
                  paste("PCC=", PCC["Limb", PCC], "| Limb obs."),
                  "Limb pred."),
  col= c("black", "grey"),
  bw.max = c(280, 150, 130, 50, 90, 50)
)
dev.off()