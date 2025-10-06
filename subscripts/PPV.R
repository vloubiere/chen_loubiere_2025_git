setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v2.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set %in% c("test", "random")]
meta[, col:= c("blue", "limegreen", "gold", 'red', 'purple', 'sienna')[.GRP], tissue]

# Plot ----
pdf("pdf/0_paper/PPV_per_tissue.pdf", width = 6.25, height = 18)
vl_par(mfrow= c(6,2))
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
  
  # Retrieve unique enhancers IDs
  dat[!grepl("^seq", ID), enh:= tstrsplit(ID, ":", keep= 1)]
  
  # Plot tile-lvl PPV
  vl_PPV(
    dat$Predictions,
    dat$score, 
    plot = T, 
    col= adjustcolor(col, .6), 
    main= paste0(tissue, " (", formatC(nrow(dat), big.mark = ","), " aug. tiles)"), 
    show.max = FALSE
  )
  
  # Add legend
  cutoff <- dat[set=="test", vl_PPV(Predictions, score, plot = T, col= col, add= T)]$predict_cutoff
  vista_total <- length(unique(na.omit(dat$enh)))
  vista_active <- length(unique(na.omit(dat[score==1, enh])))
  vista_left <- length(unique(na.omit(dat[score==1 & Predictions >= cutoff, enh])))
  vl_legend(
    "topleft",
    legend = c(
      "Test set:",
      paste(vista_total, "VISTA tiles"),
      paste(vista_active, "active enh."),
      paste(vista_left, "enh. >= cutoff")
    )
  )
  # Add legend
  vl_legend(
    legend= c(
      "Test set",
      "+ 600k rdm. seq."
    ),
    lwd= 1,
    col= c(col, adjustcolor(col, .4))
  )
  
  # Plot enhancer-level PPV
  enh <- dat[!is.na(enh), .(max.pred= max(Predictions)), .(enh, score)]
  vl_PPV(
    enh$max.pred, 
    enh$score, 
    plot = T, 
    col= "black", 
    main= paste(tissue, "(max pred. / VISTA tile)"),
    Nleft = 20
  )
  
  print(".")
}, .(tissue, col)]
dev.off()