setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Compute PCC file ----
PCC_file <- "db/PCC/PCC_ATAC_models.rds"
if(!file.exists(PCC_file)) {
  
  # Import metadata
  meta <- readRDS("Rdata/paper_metadata_v2.rds")
  meta <- meta[dataset=="accessibility" & set=="testBestDesign"]
  
  # Compute PCC
  meta[, cor := {
    # Import observed and predicted scores
    .c <- merge(fread(obs_file), fread(pred_file), by.x= "ID", by.y= "location")
    print(paste0(.GRP, "/", .NGRP))
    # Compute PCC
    cor(.c$score, .c$Predictions)
  }, .(ID, obs_file, pred_file, tissue, fold, replicate)]
  
  # Save
  saveRDS(meta, PCC_file)
}

# Import and simplify names ----
dat <- readRDS(PCC_file)
dat[, ID:= gsub("model1_bulkATAC_", "", ID)]

# Average PCC ----
av <- dcast(dat, tissue~ID, value.var = "cor", fun.aggregate = function(x) mean(x, na.rm= T))
av <- as.matrix(av, 1)
colnames(av) <- gsub("model1_bulkATAC_", "", colnames(av))
sd <- dcast(dat, tissue~ID, value.var = "cor", fun.aggregate = function(x) sd(x))
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
           pdf.file = "pdf/0_paper/heatmap_PCC_20_different_models_S1.pdf",
           pdf.cell.size = .2)

