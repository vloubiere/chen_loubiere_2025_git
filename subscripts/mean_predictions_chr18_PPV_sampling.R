setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
meta <- readRDS("Rdata/metadata_ATACSeq_model_predictions.rds")
meta <- meta[file.exists(pred_VISTA_rdm)]

# Import predicted & observed values test set ----
dat <- meta[, {
  .c <- fread(pred_VISTA_rdm)
  .c[, .(ID= location, obs= 0, pred= Predictions)]
}, .(model= ID, tissue, fold, replicate)]

# Compute mean predictions across fold/replicates ----
res <- dat[, .(obs= mean(obs), pred= mean(pred)), keyby= .(tissue, model, ID)]
res <- dcast(res, model+ID+obs~tissue, value.var = "pred")
res[, set:= "chr18_rdm"]
res[, class:= "globallyClosed"]
setcolorder(res, c("set", "class"))
saveRDS(res,
        "db/predictions/model1_VISTA_PPV_sampling_mean_predictions.rds")
