setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
meta <- readRDS("Rdata/metadata_ATACSeq_model_predictions.rds")

# Import and subset observed values chr 18 ----
obsChr18 <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/fasta/testing_dataset/whole_chr18_activity.rds")
# Remove centromeres and telomeres
heterochromatin <- fread("/groups/stark/shenzhi.chen/db/blacklisted/telomere_centromeres_mm10.bed",
                         col.names =  c("seqnames", "start", "end", "name"))
heterochromatin[, start:= start-10e3]
heterochromatin[, end:= end+10e3]
obsChr18 <- vl_intersectBed(obsChr18,
                            heterochromatin,
                            invert = T)
# Randomly sample 10,000 regions
set.seed(1)
obsChr18 <- obsChr18[sample(.N, 10000)]
obsChr18[, ID:= paste0(seqnames, ":", start, "-", end, "_", strand)]
# Melt
obsChr18 <- melt(obsChr18,
                 id.vars = "ID",
                 measure.vars = c("heart", "forebrain", "midbrain", "hindbrain", "neuralTube", "limb"),
                 variable.name = "tissue",
                 value.name = "obs")
obsChr18[, obs:= log2(obs+1)]

# Import predicted values ----
dat <- meta[, {
  .c <- fread(pred_ATAC_chr18)[, .(ID= location, pred= Predictions)]
  merge(obsChr18[.BY, on= "tissue"],
        .c,
        by = "ID")
}, .(model= ID, tissue, fold, replicate)]

# Compute delta values ----
dat[dat[tissue=="heart"], heart:= i.obs, on= c("model", "ID", "fold", "replicate")]
dat[dat[tissue=="forebrain"], forebrain:= i.obs, on= c("model", "ID", "fold", "replicate")]
dat[, obs_delta:= fifelse(tissue=="heart", obs-forebrain, obs-heart)]
dat[dat[tissue=="heart"], heart:= i.pred, on= c("model", "ID", "fold", "replicate")]
dat[dat[tissue=="forebrain"], forebrain:= i.pred, on= c("model", "ID", "fold", "replicate")]
dat[, pred_delta:= fifelse(tissue=="heart", pred-forebrain, pred-heart)]
dat$heart <- dat$forebrain <- NULL

# Melt data ---
mdat <- melt(dat,
             id.vars = c("model", "tissue", "fold", "replicate", "ID"),
             measure.vars = patterns("obs"= "^obs", "pred"= "^pred"))
mdat <- na.omit(mdat)
mdat[, variable:= c("tissue", "delta")[variable]]

# Compute PCCs ----
pcc <- mdat[, .(globPCC= .SD[, cor(obs, pred)]), .(model, tissue, fold, replicate, variable)]

# Melt for plotting ----
pl <- melt(pcc,
           id.vars = c("tissue", "model", "fold", "replicate", "variable"),
           measure.vars = "globPCC",
           variable.name = "group")
pl[, variable:= factor(variable, c("tissue", "delta"))]
setorderv(pl, c("fold", "replicate"))
pl <- pl[, {
  .(mean= mean(value), sd= sd(value), value= .(value))
}, keyby= .(variable, group, tissue, model)]

# Plot PCCs per fold and rep ----
cols <- uniqueN(pl[, tissue])
rows <- uniqueN(pl[, .(variable, group)])
ylim <- c(0, max(unlist(pl$mean)))

pdf("pdf/ATAC_PCC_fold_and_replicates_chr18.pdf", 2.5*cols, 2*rows)
vl_par(mai= c(.7, .4, .5, .2),
       mfrow= c(rows, cols))
pl[, {
  vl_barplot(height = mean,
             individual.var = value,
             bar.labels = round(mean, 2),
             bar.labels.cex = .4,
             ind.col = rep(c("tomato", "cornflowerblue", "limegreen"), each= 2),
             names.arg = gsub("^model1_bulkATAC_", "", model),
             ylim= ylim,
             ylab= paste(variable, group),
             main= tissue)
}, .(group, variable, tissue)]
dev.off()

# Compute mean predictions across fold/replicates ----
res <- mdat[, lapply(.SD, mean), .(model, tissue, variable, ID), .SDcols= c("obs", "pred")]
saveRDS(res,
        "db/predictions/model1_ATAC_mean_predictions_chr18.rds")
