setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW"]
meta <- meta[set=="test" & tissue %in% c("heart", "limb", "midbrain")]

# Predicted act at max PPV ----
meta[, predict_cutoff:= {
  
  # Import data for each set
  dat <- .SD[, {
    cmb <- .SD[, {
      if(is.na(obs_file)) {
        .c <- fread(pred_file)[, score:= 0]
        setnames(.c, "location", "ID")
        .c
      } else {
        merge(fread(obs_file), fread(pred_file), by.x= "ID", by.y= "location")
      }
    }, .(fold, replicate)]
    # SANITY CHECK -> chr18 sequences should be there 6 times and the other 2
    if(!all(unique(cmb[, .N, ID]$N) %in% c(2, 6)))
      stop("Some sequences don't have expected number of replicates")
    
    # Mean between fold/replicates
    cmb[, lapply(.SD, mean), ID, .SDcols= c("score", "Predictions")]
  }, .(set)]
  
  # Predict cutoff
  vl_PPV(
    dat$Predictions,
    dat$score,
    plot= F
  )$predict_cutoff
  
}, tissue]

# Melt ----
meta <- melt(
  meta,
  id.vars = c("tissue", "fold", "replicate", "obs_file", "predict_cutoff"),
  # measure.vars = c("pred_file", "predActFromAcc_file", "predActFromMTL_file", "predActFromRandomIni_file"),
  measure.vars = c("pred_file", "predActFromAcc_file", "predActFromRandomIni_file"),
  value.name = "pred_file"
)

# Plot ----
res <- meta[, {
  
  # Import data for each set
  cmb <- .SD[, {
    merge(fread(obs_file), fread(pred_file), by.x= "ID", by.y= "location")
  }, .(fold, replicate, variable)]
  
  # Mean between fold/replicates
  cmb <- cmb[, lapply(.SD, mean), .(ID, variable), .SDcols= c("score", "Predictions")]
  
  # Dcast
  cmb <- dcast(cmb, ID+score~variable, value.var = "Predictions")
  
  # Scale predictions from the accessibility model from 0 to 1
  cmb[, predActFromAcc_file:= (predActFromAcc_file-min(predActFromAcc_file))/diff(range(predActFromAcc_file))]
  
  # Plot PPV using randomly initialized model
  rdmInit <- vl_PPV(
    cmb$predActFromRandomIni_file, 
    cmb$score,
    return.value.at = predict_cutoff
  )
  
  # Add activity model prediction
  accOnly <- vl_PPV(
    cmb$predActFromAcc_file, 
    cmb$score,
    return.value.at = predict_cutoff
  )
  
  # Add PPV after transfer-learning 
  TL <- vl_PPV(
    cmb$pred_file, 
    cmb$score,
    return.value.at = predict_cutoff
  )
  
  # Return
  .(rdmInit, accOnly, TL)
}, .(tissue, predict_cutoff)]

# Fix labels
res[, tissue:= switch(tissue, "heart"= "Heart", "limb"= "Limb", "midbrain"= "Midbrain"), tissue]
mat <- as.matrix(res[, !"predict_cutoff"], 1)
mat <- t(mat)
mat[is.na(mat)] <- 0

# Ploting parameters
Cc <- adjustcolor(c("red", "blue", "black"), .2)

# Plot
pdf("pdf/0_paper/barplot_PPV_shuffled_models_cutoff.pdf", 3, 3)
vl_par(mai= c(1, .9, 1, 1.3))
bar <- barplot(mat,
               beside = T,
               ylab= paste("PPV at cutoff (%)"),
               col= Cc,
               xaxt= "n")
tiltAxis(x= colMeans(bar), labels= colnames(mat))
vl_legend(
  legend= c(
    "Rdm init.",
    "Acc. model (scaled)",
    "TL model"
  ),
  fill= Cc
)
dev.off()
