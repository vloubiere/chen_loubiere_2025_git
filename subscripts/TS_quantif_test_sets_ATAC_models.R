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

# Plot quantification TS groups ----
setorderv(dat,
          c("variable", "tissue"))
rows <- uniqueN(dat[, tissue])
cols <- uniqueN(dat[, variable])
pdf("pdf/quantif_TS_groups_ATAC_models_mean_predictions_test_sets.pdf",
    width= 7*cols,
    height = 3.5*rows)
vl_par(mai= c(1, .5, .5, .1),
       mfcol= c(rows, cols))
dat[grepl("specific", class, ignore.case = T), {
  # Unique vars
  models <- gsub("^model1_bulkATAC_", "", sort(unique(model)))
  classes <- unique(class)
  # pval list 
  pval <- matrix(1:(length(models)*length(classes)), ncol= 2, byrow = T)
  pval <- split(pval, seq(nrow(pval)))
  # Boxplot
  box <- vl_boxplot(pred~class+model,
                    tilt.names= T,
                    xaxt= "n",
                    col= c("cornflowerblue", "tomato"),
                    ylab= "Predicted",
                    violin= T,
                    boxwex= .3,
                    viowex= .5,
                    compute.pval= pval)
  title(main= paste(tissue, variable), line= 1.5)
  # Observed median
  sub <- .SD[model==model[1]]
  abline(h= median(sub[class=="specificClosed", pred]), col= "cornflowerblue", lty= 2)
  abline(h= median(sub[class=="tissueSpecific", pred]), col= "tomato", lty= 2)
  # axis
  vl_tilt_xaxis(x = seq(models)*length(classes)-0.5,
                labels = models)
  # Add FC
  sapply(seq(pval), function(i) {
    .c <- pval[[i]]
    diff <- box$stats[3, .c[2]]-box$stats[3,.c[1]]
    text(box$x[i],
         box$y[i],
         round(diff, 2),
         pos= 3,
         cex= 7/12,
         xpd= NA)
  })
  print(".")
}, .(variable, tissue)]
dev.off()