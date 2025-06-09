setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import metadata ----
dat <- readRDS("db/predictions/20240807_ATACSeq_model_testSet_predictions.rds")
dat <- dat[fold %in% paste0("fold0", 1:7)]

# Mean replicates ----
meanDat <- dat[, .(predicted= mean(predicted)), .(tissue, ID, class, observed)]
meanDat[!(class %in% c("globallyOpen", "globallyClosed", "specificClosed")), class:= tissue]
meanDat[, Cc:= switch(class,
                      "globallyOpen"= "red",
                      "globallyClosed"= "grey",
                      "specificClosed"= "cornflowerblue",
                      "limegreen"), class]
set.seed(1)
meanDat <- meanDat[sample(.N, .N)]

# Compute PCCs
meanDat[, Overall:= paste("PCC=", round(cor(observed, predicted), 2)), tissue]
meanDat[, GloballyOpen:=   paste("GloballyOpen=", round(cor(observed[Cc=="red"], predicted[Cc=="red"]), 2)), tissue]
meanDat[, GloballyClosed:= paste("GloballyClosed=", round(cor(observed[Cc=="grey"], predicted[Cc=="grey"]), 2)), tissue]
meanDat[, SpecificClosed:= paste("SpecificClosed=", round(cor(observed[Cc=="cornflowerblue"], predicted[Cc=="cornflowerblue"]), 2)), tissue]
meanDat[, TissueSpecific:= paste("TissueSpecific=", round(cor(observed[Cc=="limegreen"], predicted[Cc=="limegreen"]), 2)), tissue]

# Plot
pdf("pdf/ATAC_scatterplots_obs_predicted_testSet.pdf", 4, 4)
vl_par()
meanDat[, {
  plot(observed,
       predicted,
       col= adjustcolor(Cc, .4),
       pch= 16,
       cex= .2,
       main= tissue,
       xlab= "Observed (log2)",
       ylab= "Predicted (log2)")
  legend("bottomright",
         legend= c(Overall, GloballyOpen, GloballyClosed, SpecificClosed, TissueSpecific),
         text.col= c("black", "red", "grey", "cornflowerblue", "limegreen"),
         bty= "n",
         cex= 7/12)
  .SD
}, .(tissue, Overall, GloballyOpen, GloballyClosed, SpecificClosed, TissueSpecific)]
dev.off()