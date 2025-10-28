setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & tissue %in% c("heart", "limb", "midbrain")]
meta <- meta[set %in% c("test", "random", "NegGenomicRegions")]

# Predicted act at max PPV ----
meta[, c("max_PPV", "predict_cutoff"):= {
  
  # Import data for each set
  dat <- .SD[set=="test", {
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
  )[c("PPV_at_cutoff", "predict_cutoff")]
  
}, tissue]

# Compute FALSE positive ----
FP <- meta[, {
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
  
  # Check percentage of FALSE positive
  dat[score==0, .(N= .N, perc= sum(Predictions<predict_cutoff)/.N*100), set]
}, .(tissue, predict_cutoff, max_PPV)]

# Fix labels ----
FP[, `Set of regions`:= paste0(sprintf("%.2f", perc), "%")]
FP[, clean.tissue:= switch(tissue, "heart"= "Heart", "limb"= "Limb", "midbrain"= "Midbrain"), tissue]
FP[, `Set of regions`:= switch(set, "test"= "VISTA inactive seq.", "NegGenomicRegions"= "Closed genomic regions", "random"= "Random sequences"), set]
FP[, clean.tissue:= factor(clean.tissue, c("Heart", "Limb", "Midbrain"))]
FP[, `Set of regions`:= factor(`Set of regions`, c("VISTA inactive seq.", "Closed genomic regions", "Random sequences"))]
setorderv(FP, "`Set of regions`")
FP[, clean.N:= ifelse(set=="test", "~26,000", formatC(N, big.mark = ","))]
FP[, `Set of regions`:= paste0(`Set of regions`, " (n= ", clean.N, ")")]
FP[, `Set of regions`:= factor(`Set of regions`, unique(`Set of regions`))]
# Cast ----
FP[, clean.perc:= paste0(sprintf("%.2f", perc), "%")]
perc <- dcast(FP, `Set of regions`~clean.tissue, value.var= "clean.perc")
cutoff <- transpose(unique(FP[, .(clean.tissue, "Max. PPV (%)"= round(max_PPV, 1), "Act. cutoff"= round(predict_cutoff, 3))]), make.names = T, keep.names = "")

# Plot ----
pdf("pdf/0_paper/table_rejected_inactive_seq.pdf")
plotTable(cutoff, wrap.text = 23)
plotTable(perc, wrap.text = 23)
title(main= paste("Control regions with a predited activity score < cutoff"), line= -7, font.main= 1, cex.main= 1)
dev.off()