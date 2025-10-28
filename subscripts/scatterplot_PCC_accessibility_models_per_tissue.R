setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="accessibility" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set=="test"]

# Plotting parameters ----
Cc <- c("grey70", "grey20")

# Scatter plot ----
pdf("pdf/0_paper/scatterplot_PCC_accessibility_models_per_tissue.pdf", 3*3, 3*3)
vl_par(mfrow= c(3,3))
meta[, {
  
  # Import observed and predicted values
  .c <- .SD[, {
    merge(fread(obs_file), fread(pred_file), by.x= "ID", by.y= "location")
  }, .(obs_file, pred_file, replicate, fold)]
  .c[, c("ID", "label"):= tstrsplit(ID, "__")]
  .c[, label:= factor(label, c("closed", "open"))]
  
  # SANITY CHECK! ----
  # Sequences on chromosome 18 should have 6 predicted values (2 rep, 3 folds)
  # Sequences on other chromosomes should have 2 (2 replicates, no overlap betwee folds)
  .c[, Nrep:= .N, ID]
  if(any(!.c$Nrep %in% c(2, 6))) 
    warning("Some sequences have more than 6 replicates. This should not be the case and should be fixed!", call. = F)
  
  # Compute mean predicted value per sequence
  .c <- .c[, lapply(.SD, mean), .(class, ID, label), .SDcols= c("score", "Predictions")]
  
  # Add colors
  .c[, class:= factor(class, c("globallyClosed","globallyOpen","specificClosed","tissueSpecific"))]
  .c[, col:= adjustcolor(Cc[label], .5)]
  
  # Scatterplot random subset
  set.seed(1)
  .c[sample(.N, 5000), {
    vlite::rasterScatterplot(score,
                             Predictions,
                             col= col,
                             cex = .5,
                             xlab= "Observed",
                             ylab= "Predicted",
                             main= tissue)
    
  }]
  
  # Add PCC
  addPcc(cor(.c$score, .c$Predictions))
  
  # Add Densities
  setorderv(.c, "label")
  vlite::addDensity(.c$score,
                    .c$Predictions,
                    col= .c$col)
  
  # Legend
  vl_legend(legend = levels(.c$label),
            fill= Cc,
            y.adj = 2)
  
  print(".")
}, tissue]
dev.off()