setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
dat <- fread("db/predictions/model1/VISTA_mean_predictions.txt")
dat[, added:= grepl("^chr18", ID)] # Randomly sampled sequences from chromosome 18

# Compute lines and columns ----
dat[, model:= factor(model, fread("Rdata/model1_ordering.txt", header= F)[[1]])]
setorderv(dat,
          c("tissue", "model"))
dat[, line:= .GRP, tissue]
dat[, endLine:= line==.NGRP, tissue]
dat[, column:= .GRP, model]
dat[, endColumn:= line==.NGRP, model]
ncol <-  uniqueN(dat[, .(model)])
nrow <-  uniqueN(dat[, .(tissue)])

# Plot PPV
pdf("pdf/model1_VISTA_PPV_curves.pdf",
    width = 1.05*ncol,
    height = 1.05*nrow)
vl_par(mfrow= c(nrow, ncol),
       mai= c(.3, .3, 0, 0),
       omi= c(0.5, 1, .25, 0),
       font.main= 1,
       cex.main= .9,
       xpd= NA)
# dat[tissue==tissue[1] & model %in% unique(model)[1:2], {
# dat[model==model[1], {
dat[, {
  vl_PPV(predicted = pred[!added],
         label = obs[!added],
         plot = T,
         xaxt= "n",
         xlab= NA,
         ylab= NA)
  vl_PPV(predicted = pred,
         label = obs,
         plot = T,
         add = T,
         col= "red",
         pos.max= 2)
  axis(1, padj= -1)
  if(line==1)
    title(main= gsub("model1_bulkATAC_", "", model[1]),
          xpd= NA,
          line= 1)
  if(line==1 & column==1)
    legend(grconvertX(0, "ndc", "user"),
           grconvertY(1, "ndc", "user"),
           c("Test set", "Aug. test set"),
           lty= 1,
           col= c("black", "red"),
           xpd= NA,
           bty= "n",
           cex= .6)
  if(column==1)
  {
    title(ylab= "Pos. pred. value (%)",
          xpd= NA,
          line= 1)
    text(par("usr")[1],
         mean(par("usr")[3:4]),
         tissue,
         pos= 2,
         offset= 3.5,
         xpd= NA)
  }
  if(endLine)
    title(xlab= "Prediction score",
          xpd= NA,
          line= 0.5)
}, .(tissue, model, line, endLine, column, endColumn)]
dev.off()