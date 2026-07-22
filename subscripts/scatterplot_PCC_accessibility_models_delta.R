setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="accessibility" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set=="test"]
meta <- meta[tissue!="CNS"]

# Control tissue to compare
meta[, tissue_ctl:= fcase(
  tissue=="limb", "midbrain",
  tissue=="heart", "midbrain",
  default = "heart"
)]
ctl_tiss <- meta[tissue %in% c("midbrain", "heart"), !"tissue_ctl"]

# Merge
meta <- merge(
  x = meta, 
  y = ctl_tiss, 
  by.x= c("tissue_ctl", "fold", "replicate"), 
  by.y= c("tissue", "fold", "replicate"), 
  suffixes= c("", "_ctl")
)

# Plotting parameters ----
Cc <- c("tomato", "limegreen", "cornflowerblue", "grey")

# Scatter plot ----
pdf("pdf/0_paper/scatterplot_PCC_accessibility_models_delta.pdf", 3*3, 3*3)
vl_par(mfrow= c(3,3))
meta[, {
  
  # Import observed and predicted values
  .c <- .SD[, {
    # Tissue model
    tiss <- merge(fread(obs_file), fread(pred_file), by.x= "ID", by.y= "location")
    # Control tissue
    ctl <- merge(fread(obs_file_ctl), fread(pred_file_ctl), by.x= "ID", by.y= "location")
    # Split activity label from ID before merging
    tiss[, c("ID", "label"):= tstrsplit(ID, "__")]
    ctl[, c("ID", "label"):= tstrsplit(ID, "__")]
    cmb <- merge(tiss, ctl, by= "ID", suffixes= c("", "_ctl"))
    # Compute Deltas
    cmb[,
        .(
          ID,
          label,
          label_ctl,
          score_delta= score-score_ctl,
          pred_delta= Predictions-Predictions_ctl
        )
    ]
  }, .(obs_file, pred_file, obs_file_ctl, pred_file_ctl, replicate, fold)]
  
  # SANITY CHECK! ----
  # Sequences on chromosome 18 should have 6 predicted values (2 rep, 3 folds)
  # Sequences on other chromosomes should have 2 (2 replicates, no overlap betwee folds)
  .c[, Nrep:= .N, ID]
  if(any(!.c$Nrep %in% c(2, 6))) 
    warning("Some sequences have more than 6 replicates. This should not be the case and should be fixed!", call. = F)
  
  # Compute class ----
  lvl <- c("Shared", paste0(tissue, "-specific"), paste0(tissue_ctl, "-specific"), "Closed")
  .c[, class:= fcase(label=="open" & label_ctl=="open", lvl[1],
                     label=="open" & label_ctl=="closed", lvl[2],
                     label=="closed" & label_ctl=="open", lvl[3],
                     default = lvl[4])]
  .c[, class:= factor(class, lvl)]
  
  # Compute mean predicted value per sequence
  .c <- .c[, lapply(.SD, mean), .(class, ID), .SDcols= c("score_delta", "pred_delta")]
  
  # Add colors
  .c[, col:= adjustcolor(Cc[class], .5)]
  
  # Scatterplot random subset
  Ntot <- nrow(.c)
  set.seed(1)
  .c[sample(.N, 5000), {
    vlite::rasterScatterplot(
      score_delta,
      pred_delta,
      col= col,
      cex = .5,
      xlab= "Delta observed",
      ylab= "Delta predicted",
      main= paste0(tissue, "-", tissue_ctl, "\n(n= ", Ntot, ")")
    )
  }]
  
  # Add PCC
  addPcc(cor(.c$score_delta, .c$pred_delta))
  
  # Add Densities
  setorderv(.c, "class")
  vlite::addDensity(.c$score_delta,
                    .c$pred_delta,
                    col= .c$col)
  
  # Legend
  vl_legend(legend = levels(.c$class),
            fill= Cc,
            y.adj = 2)
  print(".")
  
}, .(tissue, tissue_ctl)]
dev.off()