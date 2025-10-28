setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW"]
meta <- meta[set=="test" & tissue!="CNS"]
meta <- melt(
  meta,
  id.vars = c("tissue", "fold", "replicate","obs_file" ),
  # measure.vars = c("pred_file", "predActFromAcc_file", "predActFromMTL_file", "predActFromRandomIni_file"),
  measure.vars = c("pred_file", "predActFromAcc_file", "predActFromRandomIni_file"),
  value.name = "pred_file"
)

# Plot ----
pdf("pdf/0_paper/PPV_compare_control_shuffled_models.pdf", width = 3.6*3, height = 9)
vl_par(mfrow= c(3,3), mai= c(.9, .7, .9, 1.6))
meta[, {
  
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
  ppv <- vl_PPV(
    cmb$predActFromRandomIni_file, 
    cmb$score,
    plot = T, 
    col= adjustcolor("red", .6), 
    main= paste0(tissue, " (", formatC(nrow(cmb), big.mark = ","), " aug. tiles)"), 
    show.max = FALSE,
    show.max.value = FALSE,
    show.pred.cutoff = F,
    xlim= c(0, 1),
    ylim= c(0, 100)
  )
  
  # Add activity model prediction
  ppv1 <- vl_PPV(
    cmb$predActFromAcc_file, 
    cmb$score,
    add = T, 
    col= adjustcolor("blue", .6), 
    show.max = FALSE,
    show.max.value = FALSE,
    show.pred.cutoff = F
  )
  
  # # Add swapped tissue
  # ppv2 <- vl_PPV(
  #   cmb$predActFromMTL_file, 
  #   cmb$score,
  #   add = T, 
  #   col= adjustcolor("grey30", .6), 
  #   show.max = TRUE,
  #   show.max.value = FALSE,
  #   show.pred.cutoff = F
  # )
  
  # Add PPV after transfer-learning 
  ppv3 <- vl_PPV(
    cmb$pred_file, 
    cmb$score,
    add = T, 
    col= adjustcolor("black", .6), 
    show.max = FALSE,
    show.max.value = FALSE,
    show.pred.cutoff = F
  )
  
  # Add legend
  vl_legend(
    legend= c(
      "TL model",
      "Acc. model (scaled)",
      # paste0("Swap tissue labels (max= ", round(ppv2, 1), "%)"),
      "Rdm init."
      ),
    lwd= 1,
    col= adjustcolor(c("black", "blue", "red"), 6)
    )
  print(".")
}, .(tissue)]
dev.off()