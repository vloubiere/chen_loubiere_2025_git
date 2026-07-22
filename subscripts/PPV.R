setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
# devtools::load_all("/groups/stark/vloubiere/vlite/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & tissue != "CNS"]
meta <- meta[set %in% c("test", "random", "NegGenomicRegions")]

# Plot ----
# pdf("pdf/0_paper/PPV_per_tissue.pdf", width = 6.25, height = 15)
# vl_par(mfrow= c(6,2))
meta[, {
  
  # Import data for each set ----
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
  
  # Retrieve unique enhancers IDs ----
  dat[!grepl("^seq", ID), enh:= tstrsplit(ID, ":", keep= 1)]
  dat[ID=="chr18", enh:= NA] # Should be set to NA (see enhancer-lvl PPV)
  
  # Plot tile-lvl PPV ----
  # Random sequences + test set
  dat[set %in% c("test", "random"), {
    vl_PPV(
      Predictions,
      score, 
      plot = T, 
      col= adjustcolor("red", .8),
      show.max = FALSE,
      show.pred.cutoff = F
    )
  }]
  # Randomly sampled genomic seq + test set
  dat[set %in% c("test", "NegGenomicRegions"), {
    vl_PPV(
      Predictions,
      score, 
      plot = T, 
      col= adjustcolor("blue", .8),
      show.max = FALSE,
      show.pred.cutoff = F,
      add= T
    )
  }]
  # Test set only
  cutoff <- dat[set=="test", {
    vl_PPV(
      Predictions,
      score, 
      plot = T, 
      col= adjustcolor("black", .8),
      add= T
    )
  }]
   
  bg <- dat[set=="test", sum(score)/.N]*100
  print(tissue)
  print(bg)
  print(cutoff$PPV_at_cutoff/bg)
}, tissue]
# dev.off()