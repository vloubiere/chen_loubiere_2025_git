setwd("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/clean_version/")

# Compute PCC ----
PCC_file <- "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/PCC/PCC_ATAC_models.rds"
if(!file.exists(PCC_file)) {
  
  # Import metadata
  meta <- readRDS("Rdata/data_preprocess/3.metadata_predicted_accessibility.rds")
  meta[, check:= file.exists(pred_access_test)]
  
  # Only keep existing files
  meta <- meta[(check)]
  
  # COmpute PCC
  meta[, cor := {
    obs <- fread(obs_ATAC_test)$score
    pred <- fread(pred_access_test)$Predictions
    print(.GRP)
    cor(obs, pred)
  }, .(tissue, fold, replicate, obs_ATAC_test, pred_access_test)]
  
  # Save
  saveRDS(meta,
          PCC_file)
}

# Average PCC ----
setorderv(meta, "weight")
av <- dcast(meta, tissue~ID, value.var = "cor", fun.aggregate = function(x) mean(x, na.rm= T))
av <- as.matrix(av, 1)
colnames(av) <- gsub("model1_bulkATAC_", "", colnames(av))
sd <- dcast(meta, tissue~ID, value.var = "cor", fun.aggregate = function(x) sd(x))
sd <- as.matrix(sd, 1)
sd[!is.na(sd)] <- paste0("±",round(sd[!is.na(sd)], 3))

# Heatmap ----
vl_par(mai= c(2, 1.5, .5, 1.5))
vl_heatmap(av,
           cluster.cols = T,
           col= c("blue", "white", "red"),
           legend.title = "Mean PCC",
           show.numbers = sd,
           numbers.cex = .35,
           pdf.file = "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/pdf/figures/PCC_different_models.pdf",
           pdf.cell.size = .2)

