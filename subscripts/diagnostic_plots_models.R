setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import metadata ----
dat <- readRDS("db/predictions/20240807_ATACSeq_model_testSet_predictions.rds")
dat <- dat[fold %in% paste0("fold0", 1:7)]

# Plots ----
pdf("pdf/test.pdf", width = 12, height = 4.5)
x <- c(1,2,3,4,5,5,6,7,8,9,10,11,12,13,14,15,15,15)
mat <- matrix(x, nrow= 3)
layout(mat, widths = c(1/3, 1, 1/3, 1/3, 1/3, 1))
par(omi= c(0,0,.25,0),
    cex= 1,
    cex.lab= 8/12,
    cex.axis= 7/12,
    las= 1,
    tcl= -0.1,
    mgp= c(.45, 0.2, 0),
    bty= "n",
    font.main= 1,
    pty= "s")
# meta[tissue=="Forebrain", {
meta[, {
  # Import data ----
  .c <- .SD[, fread(ATAC_predictions_testSet), .(fold, rep)]
  # .c <- .SD[, fread(ATAC_predictions_testSet, nrows= 1000), .(fold, rep)]
  .c[fread(ATAC_observed_testSet), obs:= i.score, on= "location==ID"]
  .c[, class:= tstrsplit(location, "__", keep= 2)]
  classes <- unique(.c$class)
  classesP <- seq(plyr::round_any(length(classes), 3, f= ceiling))+4
  
  # Correlations between replicates ----
  reps <- CJ(seq(unique(.c$rep)), seq(unique(.c$rep)))
  reps <- reps[V2>V1]
  reps <- reps[, lapply(.SD, function(x) paste0("rep", x))]
  reps[, {
    pcc <- merge(.c[.BY, on= "rep==V1"],
                 .c[.BY, on= "rep==V2"],
                 by= "location")
    par(mai= c(.3, .3, .3, .2))
    pcc[, {
      pl(Predictions.x,
         Predictions.y,
         xlab= paste("Predictions", V1[1]),
         ylab= paste("Predictions", V2[1]))
    }]
  }, .(V1, V2)]
  while(par("mfg")[1]<3)
    plot.new()
  
  # PCC per fold ----
  folds <- .c[, .(pcc= cor(Predictions, obs)), keyby= c("fold", "rep")]
  folds <- folds[, .(pcc= .(pcc), mean= mean(pcc), sd= sd(pcc)), fold]
  folds[, {
    par(mai= c(.5, .75, .5, .5),
        las= 2,
        pty= "m")
    bar <- barplot(mean,
                   names.arg = fold,
                   ylim= c(0, max(unlist(pcc))),
                   ylab= "pcc",
                   border= NA)
    lapply(seq(pcc), function(i) {
      points(jitter(rep(bar[i], length(pcc[[i]])), amount = .25),
             pcc[[i]],
             xpd= NA,
             cex= .5,
             col= "grey30")
    })
  }]
  
  # Mean predictions ----
  agg <- .c[, .(Predictions= mean(Predictions)), .(location, obs, class)]
  agg[, {
    par(mai= c(.5, .5, .5, .5),
        cex.main= 10/12,
        las= 1,
        pty= "s")
    pl(obs,
       Predictions,
       xlab= "Observed ATAC-Seq enrichment (log2)",
       ylab= "Average prediction (log2)")
    title(main= "Average between replicates")
  }]
  xlim <- par("usr")[1:2]
  ylim <- par("usr")[3:4]
  
  # Predictions for each class ----
  agg[, {
    par(mai= c(.3, .3, .3, .2),
        cex.main= 8/12)
    pl(obs,
       Predictions,
       xlab= "Observed (log2)",
       ylab= "Predicted (log2)",
       xlim= xlim,
       ylim= ylim)
    title(main= paste0(class, "\nn= ", formatC(.N, big.mark = ",")))
    .SD
  }, class]
  while(par("mfg")[1]*par("mfg")[2]<14)
    plot.new()
  
  # Import predictions Chr11 ----
  # .c <- .SD[fold=="fold04", fread(ATAC_predictions_testChr, nrows= 1000), .(fold, rep)]
  .c <- .SD[fold=="fold04", fread(ATAC_predictions_testChr), .(fold, rep)]
  .c <- .c[, .(Predictions= mean(Predictions)), location]
  .c[fread(ATAC_observed_testChr), obs:= log2(i.score+1), on= "location==ID"]
  
  # Plot predicted vs observed ----
  if(nrow(.c))
    .c[, {
      par(mai= c(.5, .5, .5, .5),
          cex.main= 10/12,
          las= 1,
          pty= "s")
      pl(obs,
         Predictions,
         xlab= "Observed ATAC-Seq enrichment (log2)",
         ylab= "Average prediction (log2)")
      title(main= "Chr11 (Fold 04 rep. average)")
    }]else
      plot.new()
  mtext(paste(tissue, "ATAC-Seq prediction"), outer = T, line= .25)
  
  print(paste(tissue, "DONE!"))
}, .(tissue, ATAC_observed_testSet, ATAC_observed_testChr)]
dev.off()