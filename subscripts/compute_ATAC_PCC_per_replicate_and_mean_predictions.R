setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
meta <- readRDS("Rdata/metadata_ATACSeq_model_predictions.rds")
setnames(meta, "ID", "model")
setorderv(meta, c("augmentation", "balancing", "weight"))
meta[, model:= factor(model, unique(model))]

# Import ATAC predicted & observed values test set ----
ATACtest <- meta[, {
  .c <- fread(obs_ATAC_test)
  .c[, ID:= gsub("(.*)__.*", "\\1", ID)]
  setnames(.c, "score", "obs")
  .SD[, {
    .c[, pred:= unlist(fread(pred_ATAC_test, sel= "Predictions"))]
  }, .(model, tissue, fold, replicate)]
}, obs_ATAC_test]

# Import and subset observed values chr 18 ----
obsChr18 <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/fasta/testing_dataset/whole_chr18_activity.rds")
obsChr18[, line:= .I]
# Sample 20,000 regions in the middle of the chromosome
Nsel <- 100000
obsChr18 <- obsChr18[start>40E6]
obsChr18 <- obsChr18[1:Nsel][seq(.N) %% 5==1]

# Import ATAC predicted values chrom 18 ----
ATACchr18 <- meta[, {
  # Import predictions
  .c <- fread(pred_ATAC_chr18,
              skip= obsChr18$line[1],
              nrows = Nsel)
  .c <- .c[seq(.N) %% 5==1]
  setnames(.c,
           c("ID", "pred"))
  # Add observed
  .c[, obs:= log2(obsChr18[[tissue]]+1)]
}, .(model, tissue, fold, replicate)]
ATACchr18[, class:= "chr18"]

# Combine ATAC data ----
ATAC <- rbindlist(list(test= ATACtest[, .(tissue, model, fold, replicate, class, ID, obs, pred)],
                       chr18= ATACchr18[, .(tissue, model, fold, replicate, class, ID, obs, pred)]),
                  idcol = "set")
setorderv(ATAC,
          c("set", "tissue", "model", "fold", "replicate"))

# Compute deltas ----
ATAC[ATAC[tissue=="heart"], obs_delta:= fifelse(tissue %in% c("heart", "limb"), NA_real_, obs-i.obs), on= c("model", "ID", "fold", "replicate")]
ATAC[ATAC[tissue=="hindbrain"], obs_delta:= fifelse(tissue %in% c("heart", "limb"), obs-i.obs, obs_delta), on= c("model", "ID", "fold", "replicate")]
ATAC[ATAC[tissue=="heart"], pred_delta:= fifelse(tissue %in% c("heart", "limb"), NA_real_, pred-i.pred), on= c("model", "ID", "fold", "replicate")]
ATAC[ATAC[tissue=="hindbrain"], pred_delta:= fifelse(tissue %in% c("heart", "limb"), pred-i.pred, pred_delta), on= c("model", "ID", "fold", "replicate")]

# Compute PCCs per replicate and save ----
PCC <- ATAC[, c(
  .(glob_PCC_tissue= cor(obs, pred, use = "complete.obs"),
    glob_PCC_delta= cor(obs_delta, pred_delta, use = "complete.obs")),
  .SD[grepl("specific", class, ignore.case = TRUE), 
    .(ts_PCC_tissue= ifelse(.N==0, NA_real_, cor(obs, pred, use = "complete.obs")),
      ts_PCC_delta= ifelse(.N==0, NA_real_, cor(obs_delta, pred_delta, use = "complete.obs")))]
), .(set, tissue, model, fold, replicate)]
PCC <- melt(PCC, measure.vars = patterns("PCC"))
PCC <- na.omit(PCC)
saveRDS(PCC,
        "db/statistics/model1/ATAC_PCC_replicates.rds")

# Compute mean predictions and save ----
mean <- ATAC[, lapply(.SD, mean, na.rm= T), .(set, tissue, model, class, ID), .SDcols= c("obs", "pred", "obs_delta", "pred_delta")]
fwrite(mean,
       "db/predictions/model1/ATAC_mean_predictions.txt")
