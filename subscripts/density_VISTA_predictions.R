setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
dat <- fread("db/predictions/model1/VISTA_mean_predictions.txt")
dat[, peakID:= gsub("(^.*):.*", "\\1", ID)]
dat <- dat[!grepl("^chr18", ID)] # Remove randomly sampled sequences chr. 18 (only keep selected regions)

# Add coordinates to compute ATAC-Seq overlaps ----
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
dat[vista, c("seqnames", "start", "end"):= .(seqnames, start, end), on= "peakID"]

# Add ATAC-Seq overlaps ----
peaks <- readRDS("db/peaks/bulkENCODE_confident_ATAC_peaks.rds")
peaks <- melt(peaks,
              id.vars = c("seqnames", "start", "end"),
              measure.vars = c("heart", "limb", "forebrain", "midbrain", "hindbrain", "neuralTube"),
              variable.name = "tissue")
peaks <- peaks[value==1]
dat[is.na(seqnames), ATAC:= FALSE]
dat[!is.na(seqnames), ATAC:= vl_covBed(.SD, peaks[.BY, on="tissue"])>0, tissue]

# Add groups ----
dat[, group:= fcase(ATAC & obs==1, "Open enh.",
                    !ATAC & obs==1, "Closed enh.",
                    ATAC & obs==0, "Open inact.",
                    !ATAC & obs==0, "Closed inact.")]
dat[, group:= factor(group, c("Open enh.", "Closed enh.", "Open inact.","Closed inact."))]
dat[, col:= c("tomato", "limegreen", "cornflowerblue", "grey")[group]]
dat[, model:= factor(model, fread("Rdata/model1_ordering.txt", header= F)[[1]])]
setorderv(dat,
          c("tissue", "model", "group"))

# Compute lines and columns ----
dat[, line:= .GRP, tissue]
dat[, endLine:= line==.NGRP, tissue]
dat[, column:= .GRP, model]
dat[, endColumn:= line==.NGRP, model]
ncol <-  uniqueN(dat[, .(model)])
nrow <-  uniqueN(dat[, .(tissue)])

# Plot ----
pdf("pdf/model1_VISTA_predictions_density_plots.pdf",
    width = 1*ncol,
    height = .6*nrow)
vl_par(mfrow= c(nrow, ncol),
       mai= c(.1, .1, .05, .05),
       omi= c(0.5, 1, .25, 0),
       font.main= 1,
       cex.main= .7)
dat[, {
  # if(.GRP==1)
  if(T)
  {
    # Compute densities
    lim <- range(pred, na.rm = T)
    dens <- .SD[, {
      .d <- density(pred, from= lim[1], to= lim[2])
      .(x= .(.d$x), y= .(.d$y))
    }, col]
    # Initiate plot
    plot(c(0, 1),
         c(0, max(unlist(dens$y))),
         type= "n",
         xaxt= "n")
    axis(1, pad= -1.25)
    # Plot densities
    dens[, lines(x[[1]], y[[1]], col= adjustcolor(col[1], .5), lwd= .75), col]
    # Legends
    segments(.5, 0, .5, par("usr")[4], lty= "11")
    if(line==1)
      title(main= gsub("model1_bulkATAC_", "", model[1]),
            xpd= NA,
            line= 1.5)
    if(line==1 & column==1)
      legend(grconvertX(0, "ndc", "user"),
             grconvertY(1, "ndc", "user"),
             unique(group),
             lty= 1,
             col= unique(col),
             xpd= NA,
             bty= "n",
             cex= .6)
    if(column==1)
    {
      title(ylab= "Density",
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
      title(xlab= "Predicted",
            xpd= NA,
            line= 0.5)
  }
}, .(tissue, model, line, endLine, column, endColumn)]
dev.off()