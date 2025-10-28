setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set == "test"]
meta <- meta[tissue %in% c("limb", "heart", "forebrain", "midbrain", "hindbrain", "neuralTube")]
meta <- meta[, .(fold, replicate, tissue, obs_file, pred_file)]
meta <- merge(
  meta[, .(fold, replicate, tissue, obs_file)], # Observed per tissue
  meta[tissue=="midbrain", .(fold, replicate, pred_file)], # Predicted with midbrain model
  by= c("fold", "replicate")
)

# Compute max PPV per tissue label ----
PPV <- meta[, {
  
  # Import data for each set
  .c <- .SD[, {
    .obs <- fread(obs_file)
    .pred <- fread(pred_file)
    stopifnot(nrow(.obs)==nrow(.pred))
    # Remove activity label (wrong here!)
    .obs[, ID:= tstrsplit(ID, "__", keep= 1)]
    .pred[, location:= tstrsplit(location, "__", keep= 1)]
    merge(.obs, .pred, by.x= "ID", by.y= "location")
  }, .(fold, replicate)]
  
  # Mean between fold/replicates
  .c <- .c[, lapply(.SD, mean), .(ID), .SDcols= c("score", "Predictions")]
  
  # Plot PPV using randomly initialized model
  vl_PPV(
    .c$Predictions, 
    .c$score,
    plot = F
  )$PPV_at_cutoff
  
}, .(tissue)]

# Order ----
PPV <- PPV[c("heart", "limb", "midbrain", "forebrain", "hindbrain", "neuralTube"), on= "tissue"]

# Plot ----
pdf("pdf/0_paper/PPV_midbrain_compare_other_CNS_tissues.pdf", width = 3, height = 3)
vl_par()
vl_barplot(PPV$V1,
           names.arg = PPV$tissue,
           ylab= "Max. Midbrain model PPV (%)")
dev.off()
