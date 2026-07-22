setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
meta <- readRDS("Rdata/metadata_ATACSeq_model_predictions.rds")
meta <- meta[file.exists(obs_ATAC_test) & file.exists(pred_VISTA_test)]

# Import predicted & observed values test set ----
dat <- meta[, {
  .c <- fread(obs_VISTA_test)[, .(class, ID, obs= score)]
  .SD[, {
    .c[, pred:= fread(pred_VISTA_test)$Predictions]
  }, .(model= ID, tissue, fold, replicate)]
}, obs_VISTA_test]
dat$obs_VISTA_test <- NULL
dat[, ID:= gsub("__active$|__inactive$", "", ID)] # remove tissue-specific activity labels from IDs 
dat[, peakID:= gsub("(^.*):.*", "\\1", ID)]

# Redefine classes ----
vista <- readRDS("db/folds/bulkATAC_folds.rds")[group=="vista"]
vista[, class:= fcase(rowSums(vista[, heart:neuralTube])==1, "tissueSpecific",
                      rowSums(vista[, heart:neuralTube])==6, "globallyOpen",
                      rowSums(vista[, heart:neuralTube])>0, "active",
                      default = "globallyClosed")]
dat[vista, set:= fifelse(i.seqnames=="chr18", "chr18", "test", na = "test"), on= "peakID"]
dat[vista, class:= i.class, on= "peakID"]
dat[class=="active", class:= fifelse(obs==0, "specificClosed", "shared")]

# ROC and PR AUC
AUC <- dat[, .(rocAUC= vl_ROC_AUC(predicted = pred, label = obs),
               prAUC= vl_PR_AUC(predicted = pred, label = obs)), .(set, tissue, model, fold, replicate)]
mAUC <- melt(AUC,
             measure.vars = c("rocAUC", "prAUC"))
setorderv(mAUC,
          c("set", "variable", "tissue", "model", "fold", "replicate"))
mAUC <- mAUC[, .(mean= mean(value),
                 value= .(value)), .(set, variable, tissue, model)]

# Plot PCCs per fold and rep ----
cols <- uniqueN(mAUC[, tissue])
rows <- uniqueN(mAUC[, set])

pdf("pdf/VISTA_ROC_PR_AUC_fold_and_replicates.pdf", 2.5*cols, 2*rows)
vl_par(mai= c(.7, .4, .5, .2),
       mfrow= c(rows, cols))
mAUC[, {
  vl_barplot(mean,
             bar.labels = round(mean, 2),
             bar.labels.cex = 0.4,
             individual.var = value,
             ind.col = rep(c("tomato", "cornflowerblue", "limegreen"), each= 2),
             names.arg = gsub("^model1_bulkATAC_", "", model),
             ylim= c(0, 1),
             ylab= paste(variable, set, "set"),
             main= tissue)
}, .(set, tissue, variable)]
dev.off()

# Compute mean predictions across fold/replicates ----
res <- dat[, .(obs= mean(obs), pred= mean(pred)), keyby= .(set, tissue, model, class, ID, peakID)]
saveRDS(res,
        "db/predictions/model1_VISTA_mean_predictions.rds")
