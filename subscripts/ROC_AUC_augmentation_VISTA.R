setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
source("git_deepATAC/function/compute_AUC.R")
require(pROC)

# Import data ----
dat <- readRDS("db/predictions/20240807_VISTA_model_testSet_predictions.rds")
dat[, aug:= factor(aug, c("No", "aug1", "aug2", "aug3"))]

# Compute AUC ----
ROC <- dat[, .(ROC_AUC= vl_ROC_AUC(label= active, predicted = predicted)), .(fold, tissue, rep, aug)]

# No augment only ----
noAug <- unique(ROC[aug=="No", .(tissue, ROC_AUC, rep, fold)])
noAugAUC <- noAug[, .(mean= mean(ROC_AUC), ROC_AUC= .(ROC_AUC)), tissue]

# Plot ----
pdf("pdf/VISTA_compare_folds_reps_ROC_AUC.pdf", 4.5, 3)
vl_par(mai= c(.7, .9, .9, .2))
# ROC AUC per tissue and augmentation
vl_boxplot(ROC_AUC~aug+tissue,
           ROC,
           outline= T,
           tilt.names= T,
           ylab= "ROC AUC")
abline(h= .5, lty= 3)
# ROC AUC per fold/rep no augmentation
vl_par(mai= c(.9, 1.4, .9, 1.4))
vl_barplot(noAugAUC$mean,
           individual.var = noAugAUC$ROC_AUC,
           main= "ROC AUC per rep/fold\nno augmentation",
           names.arg = noAugAUC$tissue, 
           border= NA, 
           ind.col = adjustcolor("grey10", .5),
           ylab= "ROC AUC")
abline(h= 0.5, lty= 3)
# ROC curves not augmentation
par(mfrow= c(2,3),
    mai= c(.2,.2,0,0),
    omi= c(.2, .5, .3, .05))
dat[aug=="No", {
  plot(c(0, 1),
       c(0, 1),
       xlab= c("False positive rate"),
       ylab= c("True positive rate"),
       type= "n")
  text(par("usr")[1],
       par("usr")[4]-strheight("M"),
       tissue,
       cex= 9/12,
       pos= 4)
  .SD[, {
    vl_ROC_AUC(label= active,
               predicted = predicted,
               plot.line = T)
  }, .(fold, rep)]
  abline(0, 1, lty= 3)
}, tissue]
text(grconvertX(.5, "ndc", "user"),
     grconvertY(.5, "line", "user"),
     "False positive rate",
     pos = 3,
     xpd= NA)
text(grconvertX(3, "line", "user"),
     grconvertY(.5, "ndc", "user"),
     "True positive rate",
     pos = 3,
     xpd= NA,
     srt= 90)
mtext("No augmentation",
      outer = T)
dev.off()