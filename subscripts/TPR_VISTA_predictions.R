setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
dat <- fread("db/predictions/model1/VISTA_mean_predictions.txt")
dat[, peakID:= gsub("(^.*):.*", "\\1", ID)]
dat[, added:= grepl("^chr18", ID)] # Extra negative control sequences from chr18

# Compute lines and columns ----
dat[, model:= factor(model, fread("Rdata/model1_ordering.txt", header= F)[[1]])]
setorderv(dat,
          c("tissue", "model"))
dat[, column:= .GRP, model]
ncol <-  uniqueN(dat[, .(model)])
nrow <-  uniqueN(dat[, .(tissue)])

# Compute PPV cutoff
dat[, cutoff:= vl_PPV(pred, obs)$predict_cutoff, .(tissue, model)]

# Collapse per VISTA enhancer
coll <- dat[(!added), .(pred= max(pred)), .(model, tissue, cutoff, peakID, obs, column)]

# Plot confusion matrices ----
pdf("pdf/model1_VISTA_TPR_tile_counts.pdf",
    width = 1.25*ncol,
    height = .65*nrow*2)
vl_par(mfrow= c(nrow*2, ncol),
       mai= c(.1, .15, .2, .2),
       omi= c(0.5, 1, .25, 0),
       font.main= 1,
       cex.main= .9,
       xpd= NA,
       cex= .2)
par(cex= .4)
# Count tiles ----
dat[, {
# dat[tissue=="forebrain" &  column==1, {
  # Confusion matrix
  .t <- table(pred= pred>cutoff,
              obs= obs==1)
  .t <- formatC(.t, big.mark = ",")
  .t <- matrix(.t, ncol= 2, dimnames = dimnames(.t))
  .t <- as.data.table(.t, keep.rownames = "")
  vl_plot_table(.t)
  # Legends
  text(par("usr")[1], 1, paste("Pred. >", round(cutoff, 2)), xpd= NA, srt= 90, pos= 2, offset= 2, cex= .8)
  text(0.5, 0, "Enhancer activity label", xpd= NA, pos= 1, cex= .8)
  title(main= gsub("model1_bulkATAC_", "", model), line= .5, cex= 1.5)
  if(column==1)
    text(0, .5, paste0(tissue, "\nN tiles"), pos= 2, xpd= NA, offset= 3.5, cex= 1/0.4)
}, .(tissue, model, cutoff, column)]

# Count VISTA enhancers ----
# coll[tissue=="forebrain" &  column==1, {
coll[, {
  # Confusion matrix
  .t <- table(pred= pred>cutoff,
              obs= obs==1)
  Ntiles <- unlist(.t[2,2])
  .t <- formatC(.t, big.mark = ",")
  .t <- matrix(.t, ncol= 2, dimnames = dimnames(.t))
  .t <- as.data.table(.t, keep.rownames = "")
  vl_plot_table(.t)
  # Legends
  text(par("usr")[1], 1, paste("Pred. >", round(cutoff, 2)), xpd= NA, srt= 90, pos= 2, offset= 2, cex= .8)
  text(0.5, 0, "Enhancer activity label", xpd= NA, pos= 1, cex= .8)
  title(main= gsub("model1_bulkATAC_", "", model), line= 0.5, cex= 1.5)
  if(column==1)
    text(0, .5, paste0(tissue, "\nN enh."), pos= 2, xpd= NA, offset= 3.5, cex= 1/0.4)
}, .(tissue, model, cutoff, column)]
dev.off()