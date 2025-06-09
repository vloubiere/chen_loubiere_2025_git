setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import metadata ----
meta <- readRDS("Rdata/metadata_ATACSeq_model_predictions.rds")
meta[, model:= paste0(ID, "__", tissue, "__", fold, "__", replicate)]
meta <- melt(meta,
             id.vars = "model",
             measure.vars = patterns("^pred_"),
             variable.name = "type")
meta[, c("model", "tissue", "fold", "rep"):= tstrsplit(model, "__")]
meta[, c("type", "set"):= tstrsplit(type, "_", keep= c(2,3))]

# Import folds ----
dat <- meta[, {
  .c <- .SD[, {
    fread(value, col.names = c("ID", "pred"))
  }, model]
  dcast(.c,
        ID~model,
        value.var = "pred",
        sep= "__")
}, .(type, set, tissue, fold, rep)]
dat[, ID:= tstrsplit(ID, "__", keep= 1)]

# Compute mean ----
mean <- dat[, lapply(.SD, function(x) mean(x, na.rm= T)), .(type, set, tissue, ID), .SDcols= patterns("^model1")]
mean <- dcast(mean,
              type+set+ID~tissue,
              value.var= grep("^model1", names(dat), value = T),
              sep="__")
mean[grepl("^chr", ID), "seqnames":= tstrsplit(ID, ":", keep= 1, type.convert = T)]
mean[grepl("^chr", ID), c("start", "end"):= tstrsplit(ID, ":|-|_", keep= c(2, 3), type.convert = T)]
mean[grepl("^chr", ID), c("seqnames", "start", "end", "strand"):= tstrsplit(ID, ":|-|_", keep= 1:4, type.convert = T)]
mean[!grepl("^seq", ID), strand:= ifelse(grepl("+", fixed = T, ID), "+", "-")]
mean[, head(.SD), .(type, set)][, c(1,2,3,124,125,126, 127)]

# Collapse per tile ID ----
coll <- dcast(dat,
              type+set+ID~tissue+fold+rep,
              value.var= grep("^model1", names(dat), value = T),
              sep="__")