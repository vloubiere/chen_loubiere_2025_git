setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# List peak files (reproducible from rep1 and rep2, see paper) ----
meta <- data.table(
  file= list.files(
    path = "/groups/stark/shenzhi.chen/projects/accessibility_model_enhancer_design_17112025/db/Chip/",
    pattern = "narrowPeak.gz$",
    recursive = T,
    full.names = T
  )
)

# Parse ----
meta[, c("cdition", "tissue"):= tstrsplit(file, "/", keep = 10:11)]

# List bw files ----
# meta[, bw:= list.files(
#   path = paste0("/groups/stark/shenzhi.chen/projects/accessibility_model_enhancer_design_17112025/db/Chip/", cdition, "/", tissue, "/"),
#   pattern = ".bigWig.gz$",
#   recursive = T,
#   full.names = T
# ), .(cdition, tissue)]

# Select relevant tissues ----
meta <- meta[tissue %in% c("midbrain", "limb", "heart")]

# Import all peaks ----
peaks <- meta[, importBed(file), .(tissue, cdition)]

# Save ----
saveRDS(peaks, "Rdata/revision_ChIPseq_peaks.rds")