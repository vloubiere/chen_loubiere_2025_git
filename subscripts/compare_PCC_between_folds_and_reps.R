setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import metadata ----
dat <- readRDS("db/predictions/20240807_ATACSeq_model_testSet_predictions.rds")
dat <- dat[fold %in% paste0("fold0", 1:7)]

# Consistency between replicates/folds ----
pcc <- dat[, .(pcc= cor.test(observed, predicted)$estimate), .(tissue, fold, rep)]
pcc <- pcc[, .(mean_pcc= mean(pcc), pcc= .(pcc)), tissue]

# Plot ----
pdf("pdf/ATAC_compare_folds_reps.pdf", 3.5, 3)
vl_par()
pcc[, {
  vl_barplot(mean_pcc,
             individual.var = pcc,
             main= "PCC per rep/fold",
             bar.labels = round(mean_pcc, 2),
             names.arg = tissue, border= NA, 
             ind.col = adjustcolor("grey10", .5),
             ylab= "PCC")
  legend(par("usr")[2]-strwidth("M"),
         par("usr")[4],
         legend= "3 rep. x 10 folds\n(n=30)",
         bty= "n",
         xpd= T,
         cex= 7/12)
}]
dev.off()