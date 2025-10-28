setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set == "test"]
meta <- meta[tissue %in% c("limb", "heart", "forebrain", "midbrain", "hindbrain", "neuralTube")]
meta <- meta[, .(fold, replicate, tissue, obs_file, pred_file)]
meta <- merge(
  meta[, .(fold, replicate, tissue, obs_file)], # Observed per tissue
  meta[, .(fold, replicate, tissue, pred_file)], # Predicted with midbrain model
  by= c("fold", "replicate"),
  suffixes= c(".lab", ".pred"),
  allow.cartesian= TRUE
)

# Compute max PPV per tissue label ----
max.PPV <- meta[, {
  
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
  
  # Max predicted value per vista tile
  .c[, ID:= tstrsplit(ID, ":", keep= 1)]
  .c[, .SD[which.max(Predictions)], ID]
}, .(tissue.lab, tissue.pred)]

# Boxplots ----
# max.PPV[, {
#   vl_boxplot(Predictions~tissue.obs, .SD[score==1], tilt.names= T, main= tissue.pred)
#   print(".")
# }, tissue.pred]

# Compute median ----
max.PPV[, tissue.lab:= factor(tissue.lab, c("heart", "limb", "midbrain", "forebrain", "hindbrain", "neuralTube"))]
max.PPV[, tissue.pred:= factor(tissue.pred, c("heart", "limb", "midbrain", "forebrain", "hindbrain", "neuralTube"))]
agg <- max.PPV[score==1, median(Predictions), .(tissue.lab, tissue.pred)]

# Dcast ----
mat <- dcast(agg, tissue.lab~tissue.pred, value.var = "V1")
mat <- as.matrix(mat, 1)

vl_heatmap(mat,
           cluster.rows = F,
           show.numbers = round(mat, 2),
           numbers.cex = 0.7*5/8,
           legend.title = "Median pred. value\n(max. per tile)",
           pdf.file = "pdf/0_paper/heatmap_median_pred.pdf",
           pdf.cell.size = .2)


