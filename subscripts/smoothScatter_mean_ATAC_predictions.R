setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
dat <- fread("db/predictions/model1/ATAC_mean_predictions.txt")
dat <- melt(dat,
            measure.vars = patterns("obs"= "^obs", "pred"= "^pred"))
dat[, variable:= c("tissue", "delta")[variable]]
dat <- dat[set=="test" | variable=="tissue"]
dat[, set:= switch(set,
                   "chr18"= "Rdm chr18",
                   "test"= "Sel. regions"), set]
dat[, set:= factor(set, c("Sel. regions", "Rdm chr18"))]
dat[, model:= factor(model, fread("Rdata/model1_ordering.txt", header= F)[[1]])]
setorderv(dat,
          c("set", "variable", "tissue", "model", "class"))

# Compute PCC ----
dat <- na.omit(dat)
dat[, c("gPCC", "tsPCC"):= {
  ts <- grepl("specific", class, ignore.case = T)
   .(round(cor(obs, pred), 2),
     ifelse(any(ts), round(cor(obs[ts], pred[ts]), 2), NA_real_))
}, .(set, variable, tissue, model)]

# Subset for plotting ----
dat[, sel:= {
  uniq <- unique(.SD[, .(class, ID)])
  uniq[, prob:= nrow(uniq)/.N, class]
  set.seed(1)
  ID %in% sample(uniq$ID, 500, prob= uniq$prob)
}, .(set, tissue)]

# Compute plotting parameters ----
dat[, line:= .GRP, tissue]
dat[, endLine:= line==max(line)]
dat[, column:= .GRP, model]
dat[, endColumn:= column==max(column)]
dat[, col:= switch(class,
                   "tissueSpecific"= "tomato",
                   "specificClosed"= "cornflowerblue",
                   "globallyOpen"= "grey80",
                   "grey50"), class]
ncol <-  uniqueN(dat[, .(model)])
nrow <-  uniqueN(dat[, .(tissue)])

# Plot ----
pdf("pdf/model1_ATAC_smoothScatter_mean_predictions.pdf",
    width = 1.25*ncol,
    height = 1.25*nrow)
vl_par(mfrow= c(nrow, ncol),
       mai= c(.1, .15, .2, .2),
       omi= c(0.5, 1, .25, 0),
       font.main= 1,
       cex.main= .9)
# dat[set==set[1] & tissue %in% unique(tissue)[1:2] & variable==variable[1], {
dat[, {
  # Scatterplot
  .SD[(sel), {
    vl_densScatterplot(obs,
         pred,
         label = class,
         pch= 16,
         col= adjustcolor(col, .4),
         cex= .6,
         xaxt= "n",
         dens.lw = .35,
         plot.legend = F)
    axis(1, padj= -1.25)
    abs <- max(quantile(c(obs, pred), .005, na.rm= T))
    abe <- min(quantile(c(obs, pred), .995, na.rm= T))
    segments(abs, abs, abe, abe, lty= "11")
  }]
  # Legends
  if(line==1)
    title(main= gsub("model1_bulkATAC_", "", model[1]),
          xpd= NA,
          line= 1.5)
  if(column==1)
  {
    title(ylab= "Predicted",
          xpd= NA,
          line= 1)
    text(par("usr")[1],
         mean(par("usr")[3:4]),
         paste0(tissue, "\n", set, "\n", variable),
         pos= 2,
         offset= 3.5,
         xpd= NA)
  }
  if(endLine)
    title(xlab= paste("Observed", variable),
          xpd= NA,
          line= 0.5)
  leg <- paste("glob=", gPCC)
  if(!is.na(tsPCC))
    leg <- c(leg, paste("ts=", tsPCC))
  legend("bottomright",
         title = "PCC",
         legend= leg,
         bty= "n",
         cex= .7)
  # Return
  print(paste(.GRP, "/", .NGRP))
}, .(set, variable, tissue, model, gPCC, tsPCC, line, endLine, column, endColumn)]
dev.off()