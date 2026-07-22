setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
meta <- readRDS("Rdata/metadata_ATACSeq_model_predictions.rds")
setnames(meta, "ID", "model")

# Import VISTA predicted values and add random sequences ----
VISTA <- meta[file.exists(pred_VISTA_test), {
  # Obs
  .c <- fread(obs_VISTA_test)[, .(class, ID, obs= score)]
  .c[, ID:= gsub("(.*)__.*", "\\1", ID)]
  .SD[, {
    # Pred
    .c[, pred:= fread(pred_VISTA_test, sel= "Predictions")]
    # Add globally closed random sequences from chr18
    add <- fread(pred_VISTA_chr18)
    setnames(add, c("ID", "pred"))
    add[, c("class", "obs"):= .("globallyClosed", 0)]
    rbind(.c, add)
  }, .(model, tissue, fold, replicate)]
}, obs_VISTA_test]

# Redefine classes ----
VISTA[, Ntissues:= length(unique(tissue[obs==1])), ID]
VISTA[, class:= fcase(Ntissues==1 & obs, "tissueSpecific",
                      Ntissues>4, "globallyOpen",
                      Ntissues>0 & obs==0, "specificClosed",
                      Ntissues>0 & obs==1, "active",
                      Ntissues==0, "globallyClosed")]

# Add set ----
VISTA[, peakID:= gsub("(.*):.*", "\\1", ID)]
chr18 <- readRDS("db/peaks/vista_tiles_clean.rds")[seqnames=="chr18", peakID]
VISTA[, set:= "test"]
# VISTA[grepl("^seq", ID), set:= "rdm"] # Random seq not used
VISTA[grepl("^chr18", ID) | peakID %in% chr18, set:= "chr18"]

# Compute ROC/PR AUC and mPCC
setorderv(VISTA,
          c("set", "model", "tissue", "fold", "replicate"))
stats <- VISTA[!grepl("^chr18", ID), { # Exclude randomly sampled sequences from chr18
  # Compute ROC/PR AUC and Matthew's PCC per replicate
  c(list(rocAUC= vl_ROC_AUC(predicted = pred, label = obs),
         prAUC= vl_PR_AUC(predicted = pred, label = obs),
         mcPCC= vl_mPCC(predicted = pred, label = obs),
         PPV_at_5= sum(pred>.5 & obs==1)/sum(pred>.5 | obs==1)),
    vl_PPV(predicted = pred, label = obs))
}, .(set, tissue, model, fold, replicate)]
stats <- melt(stats,
              id.vars = c("set", "tissue", "model", "fold", "replicate"))
saveRDS(stats,
        "db/statistics/model1/VISTA_AUC_mPCC_PPV_replicates.rds")

# Compute mean predictions and save ----
mean <- VISTA[, lapply(.SD, mean), .(set, tissue, model, class, ID), .SDcols= c("obs", "pred")]
fwrite(mean,
       "db/predictions/model1/VISTA_mean_predictions.txt")
