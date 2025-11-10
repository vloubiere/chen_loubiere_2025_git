setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set == "test"]
meta <- meta[tissue %in% c("limb", "heart", "forebrain", "midbrain", "hindbrain", "neuralTube")]
meta <- meta[, .(fold, replicate, tissue, obs_file, pred_file)]
meta <- merge(
  meta[, .(fold, replicate, tissue, obs_file)], # Observed per tissue
  meta[tissue=="midbrain", .(fold, replicate, tissue, pred_file)], # Predicted with midbrain model
  by= c("fold", "replicate"),
  suffixes= c(".lab", ".pred"),
  allow.cartesian= TRUE
)

# Compute max PPV per tissue label ----
max.pred <- meta[, {
  
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

# Retrieve pan-CNS labels ----
panCNS <- readRDS("db/folds/bulkATAC_folds.rds")
panCNS <- panCNS[group=="vista"]
panCNS <- panCNS[midbrain==1 & hindbrain==1 & neuralTube==1 & forebrain==1, peakID]
max.pred[, panCNS:= ID %in% panCNS]

# Rbind
sel <- list(
  max.pred[tissue.lab=="midbrain" & ID %in% max.pred[, all(score==0), ID][(V1), ID]][, tissue.lab:= "Inactive (MB)"],
  max.pred[score==1],
  max.pred[tissue.lab=="midbrain" & (panCNS)][, tissue.lab:= "Pan-CNS"]
)
sel <- rbindlist(sel)
sel[, tissue.lab:= factor(tissue.lab, c("Inactive (MB)", "heart", "limb", "midbrain", "forebrain", "hindbrain", "neuralTube", "Pan-CNS"))]

# Boxplots ----
pdf("pdf/0_paper/boxplot_midbrain_pred_act_VISTA_enhancers.pdf", width = 3, height = 3)
vl_par()
vl_boxplot(
  Predictions~tissue.lab,
  sel,
  tilt.names= T,
  compute.pval= list(c(4,8)),
  ylab= "Max. activity score"
  )
dev.off()

# # Compute median ----
# max.PPV[, tissue.pred:= factor(tissue.pred, c("heart", "limb", "midbrain", "forebrain", "hindbrain", "neuralTube"))]
# agg <- max.PPV[score==1, median(Predictions), .(tissue.lab, tissue.pred)]
# 
# # Dcast ----
# mat <- dcast(agg, tissue.lab~tissue.pred, value.var = "V1")
# mat <- as.matrix(mat, 1)
# 
# vl_heatmap(mat,
#            cluster.rows = F,
#            show.numbers = round(mat, 2),
#            numbers.cex = 0.7*5/8,
#            legend.title = "Median pred. value\n(max. per tile)",
#            pdf.file = "pdf/0_paper/heatmap_median_pred.pdf",
#            pdf.cell.size = .2)
# 
# 
