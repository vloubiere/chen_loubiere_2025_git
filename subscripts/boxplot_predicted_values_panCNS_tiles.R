setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set == "test"]
meta <- meta[tissue=="midbrain"]

# Import obs/predicted values ----
dat <- meta[, {
  
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
  .c[, lapply(.SD, mean), .(ID), .SDcols= c("score", "Predictions")]
  
}, .(tissue)]

# Max predicted activity per enhancer ID ----
dat[, ID:= tstrsplit(ID, ":", keep= 1)]
dat <- dat[, .SD[which.max(Predictions)], ID]

# Retrieve pan-CNS labels ----
panCNS <- readRDS("db/folds/bulkATAC_folds.rds")
panCNS <- panCNS[group=="vista"]
panCNS <- panCNS[midbrain==1 & hindbrain==1 & neuralTube==1 & forebrain==1, peakID]
dat[, panCNS:= ID %in% panCNS]

# Plot ----
pdf("pdf/0_paper/boxplot_pred_act_midbrain_vs_panCNS_enhancers.pdf", width = 2.3, height = 3)
vl_par()
vl_boxplot(
  list(
    `Act. in midbrain`= dat[score==1, Predictions],
    `Pan-CNS`= dat[(panCNS), Predictions]
  ),
  ylab= "Max. predicted activity",
  tilt.names= T,
  ylim= c(0, 1),
  compute.pval= list(c(1,2))
)
dev.off()