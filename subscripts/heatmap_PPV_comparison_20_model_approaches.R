setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Compute PPV file ----
PPV_file <- "db/PPV/PPV_ACTIVITY_models.rds"
if(!file.exists(PPV_file)) {
  
  # Import metadata
  meta <- readRDS("Rdata/paper_metadata_v2.rds")
  meta <- meta[dataset=="activity" & set=="test"]
  
  # Compute PV
  dat <- meta[, {
    # Import observed and predicted scores
    .c <- .SD[, {
      merge(fread(obs_file), fread(pred_file), by.x= "ID", by.y= "location")
    }, .(obs_file, pred_file, fold, replicate)]
    print(paste0(.GRP, "/", .NGRP))
    # Compute mean signal per region
    .c <- .c[, lapply(.SD, mean), ID, .SDcols= c("score", "Predictions")]
    # Compute PPV
    .(PPV= vl_PPV(.c$Predictions, .c$score, plot= F)$PPV_at_cutoff)
  }, .(ID, tissue)]
  
  # Save
  saveRDS(dat, PPV_file)
} else
  dat <- readRDS(PPV_file)
  
# Simplify names ----
dat[, ID:= gsub("model1_bulkATAC_", "", ID)]

# Average PCC ----
av <- dcast(dat, tissue~ID, value.var = "PPV")
av <- as.matrix(av, 1)
colnames(av) <- gsub("model1_bulkATAC_", "", colnames(av))

# Heatmap ----
vl_par(mai= c(2, 1.5, .5, 1.5))
vl_heatmap(av,
           cluster.cols = T,
           col= c("blue", "white", "red"),
           legend.title = "Max PPV",
           show.numbers = round(av, 1),
           numbers.cex = .35,
           pdf.file = "pdf/0_paper/heatmap_PPV_20_different_models_S1.pdf",
           pdf.cell.size = .2)

