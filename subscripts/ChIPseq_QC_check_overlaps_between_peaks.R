setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# Import peaks ----
peaks <- readRDS("Rdata/revision_ChIPseq_peaks.rds")

# Import peaks ----
vl_par(mfrow= c(2,2))
dat <- peaks[signalValue>5]
dat[, {
  .c <- collapseBed(.SD)
  .c[, name:= paste0("peak_", .I)]
  .p <- split(.SD[, !"tissue"], tissue)
  ov <- sapply(.p, function(x) intersectBed(.c, x)[, name])
  upsetPlot(
    ov
  )
  title(main= cdition[1])
  print("")
}, cdition]
