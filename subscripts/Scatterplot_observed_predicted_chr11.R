setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import metadata ----
dat <- readRDS("db/predictions/20240807_ATACSeq_model_chr11_predictions.rds")
dat <- dat[fold %in% paste0("fold0", 1:7)]

# Mean replicates ----
meanDat <- dat[, .(predicted= mean(predicted)), .(tissue, ID, observed)]

# Plot ----
pdf("pdf/ATAC_scatterplots_obs_predicted_chr11.pdf", 3, 3)
vl_par()
meanDat[, {
  obs <- log2(observed+0.01)
  smoothScatter(obs,
                predicted,
                colramp = colorRampPalette(c("white", blues9[-c(1,2)])),
                main= paste(tissue, " chr. 11"),
                xlab= "Log2(Observed+1)",
                ylab= "Predicted (log2)")
  vl_plot_coeff("bottomright", 
                cor(obs, predicted))
  # legend("bottomright",
  #        legend= c(Overall, GloballyOpen, GloballyClosed, SpecificClosed, TissueSpecific),
  #        text.col= c("black", "red", "grey", "cornflowerblue", "limegreen"),
  #        bty= "n",
  #        cex= 7/12)
  .SD
}, tissue]
dev.off()