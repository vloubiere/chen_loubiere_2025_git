setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import data ----
dat <- readRDS("db/predictions/model1_VISTA_mean_predictions.rds")
dat$peakID <- NULL
samp <- readRDS("db/predictions/model1_VISTA_PPV_sampling_mean_predictions.rds")
setkeyv(samp, c("model"))

# Add randomly sampled regions so that 1% of each set is active ----
aug <- dat[, {
  Nsamp <- round(sum(obs)*100)-.N
  set.seed(1)
  add <- samp[.BY, on= "model"]
  add <- add[sample(.N, Nsamp, replace = Nsamp>.N), c("set", "ID", "class", "obs", tissue), with= F]
  setnames(add, tissue, "pred")
  rbind(.SD,
        add)
}, .(tissue, model)]

# Plot PPV
setorderv(aug,
          c("tissue", "model"))
rows <- uniqueN(dat[, tissue])
cols <- uniqueN(dat[, model])
pdf("pdf/VISTA_PPV_test_set_only_and_rdm.pdf",
    width = cols*2.5, 
    height = rows*2.5)
vl_par(mai= rep(.5, 4),
       mfrow= c(rows, cols))
aug[, {
  # PPV only on test set
  vl_PPV(pred[set!="chr18_rdm"], obs[set!="chr18_rdm"], plot= T)
  title(main= paste(gsub("^model1_bulkATAC_", "", model), tissue))
  
  # PPV after adding random sequences
  PPV <- vl_PPV(pred, obs, plot= T, add= T, col= "red")
  legend("topleft",
         legend= c("Test set", "Add rdm"),
         lty= 1,
         col= c("black", "red"),
         bty= "n",
         cex= .7)
  PPV
}, .(tissue, model)]
dev.off()

# Median PPV per tissue ----
PPV <- aug[, {
  # PPV after adding random sequences
  c(vl_PPV(pred, obs, plot= F), list(n= sum(obs)))
}, .(tissue, model)]
PPV[, tissue:= paste0(tissue, " (n= ", formatC(n, big.mark = ","), ")"), .(tissue, n)]

# Compute median PPV for odering ----
PPV[, median_PPV:= median(PPV_at_cutoff), model]
pl <- dcast(PPV,
            tissue~median_PPV+model,
            value.var = "PPV_at_cutoff")
setnames(pl,
         gsub("^.*bulkATAC_(.*)", "\\1", names(pl)))
ylim <- range(unlist(pl[, -1]))
pdf("pdf/VISTA_PPV_compare_models_tissues.pdf", 7, 4)
vl_par(mai= c(1.3, .9, .2, 2))
pl[, {
  if(.GRP==1)
  {
    plot(NA,
         xlim= c(0, length(.SD)),
         ylim= ylim,
         xaxt= "n",
         xlab= NA,
         ylab= "PPV at cutoff (%)")
    title(xlab= "Model", line = 4.5)
    vl_tilt_xaxis(seq(.SD), labels = gsub("^model1_bulkATAC_", "", names(.SD)))
  }
  col <- rainbow(.NGRP)[.GRP]
  lines(unlist(.SD),
        col= adjustcolor(col, .3), lwd= 2)
  points(unlist(.SD),
         pch= 16,
         col= adjustcolor(col, .5))
  x <- unlist(.SD)
  max <- which.max(x)
  closeMax <- between(x, .95*max(x), max(x), incbounds = F)
  points(which(closeMax),
         x[closeMax],
         pch= 4,
         cex= .3)
  points(max,
         x[max])
}, tissue]
legend(par("usr")[2],
       par("usr")[4],
       legend= c("Maximum PPV", "PPV>95% maxPPV"),
       pch= c(1, 4),
       bty= "n",
       xpd= T,
       cex= .7)
legend(par("usr")[2],
       par("usr")[4]-strheight("M")*3,
       legend= pl$tissue,
       col= rainbow(uniqueN(pl$tissue)),
       lty= 1,
       bty= "n",
       xpd= T,
       cex= .7)
dev.off()
