setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/predictions/model1_ATAC_mean_predictions_test_set.rds")

# Compute probs for random sampling ----
dat[, tot:= .N, .(tissue, model, variable)]
dat[, prob:= tot/.N, .(tissue, model, variable, class)]

# Add colors ----
dat[, col:= switch(class, 
                   "globallyOpen"= "grey80",
                   "globallyClosed"= "grey50",
                   "specificClosed"= "cornflowerblue",
                   "tissueSpecific"= "tomato"), class]

# Plot smoothScatter plots ----
setorderv(dat,
          c("variable", "tissue", "model"))
rows <- uniqueN(dat[, tissue])
cols <- uniqueN(dat[, model])
pdf("pdf/smoothScatter_ATAC_models_mean_predictions_test_sets.pdf", width = 3*cols, height = 3*rows)
vl_par(mfrow= c(rows, cols))
dat[, {
  if(T)
  {
    # Compute PCC
    globPCC <- cor(obs, pred)
    globPCC <- round(globPCC, 2)
    tsPCC <- cor(obs[grepl("specific", class, ignore.case = T)],
                   pred[grepl("specific", class, ignore.case = T)])
    tsPCC <- round(tsPCC, 2)
    
    # Scatterplot
    set.seed(1)
    sub <- .SD[sample(.N, 2000, prob= prob)]
    sub[, {
      vl_rasterScatterplot(obs,
                           pred,
                           xlab= paste("Observed", variable),
                           ylab= paste("Predicted", variable),
                           col= adjustcolor(col, .3),
                           cex= .5,
                           pch= 16)
      abline(0, 1, lty= 3)
      legend("bottomright",
             bty= "n",
             title= "PCC",
             legend= c(paste("Global", globPCC),
                       paste("Tissue-specific", tsPCC)),
             cex= 7/12)
    }]
    title(main= paste(tissue, gsub("^model1_bulkATAC_", "", model)),
          line= 3.5)
    
    # Density
    lineW <- 1.5
    .SD[, {
      # Observed
      dens <- density(obs, from= par("usr")[1], to= par("usr")[2])
      x <- dens$x
      y <- dens$y
      y <- y/max(y)*diff(grconvertY(c(0,lineW), "line", "user"))+par("usr")[4]
      lines(x, y, col= col[1], xpd= NA)
      # Predicted
      dens <- density(pred, from= par("usr")[3], to= par("usr")[4])
      x <- dens$x
      y <- dens$y
      y <- y/max(y)*diff(grconvertX(c(0,lineW), "line", "user"))+par("usr")[2]
      lines(y, x, col= col[1], xpd= NA)
    }, col]
  }
}, .(variable, tissue, model)]
dev.off()