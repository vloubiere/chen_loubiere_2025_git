setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/predictions/model1_ATAC_mean_predictions_chr18.rds")
dat[, col:= "lightgrey"]

# Plot smoothScatter plots ----
setorderv(dat,
          c("variable", "tissue", "model"))
rows <- uniqueN(dat[, tissue])
cols <- uniqueN(dat[, model])
pdf("pdf/smoothScatter_ATAC_models_mean_predictions_chr18.pdf", width = 3*cols, height = 3*rows)
vl_par(mfrow= c(rows, cols))
dat[, {
  # if(.GRP==1)
  if(T)
  {
    # Compute PCC
    globPCC <- cor(obs, pred)
    globPCC <- round(globPCC, 2)
    
    # Scatterplot
    set.seed(1)
    sub <- .SD[sample(.N, 2000)]
    .SD[, {
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
             legend= paste("Global", globPCC),
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