setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import edit
edit <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/Leididi_design/heart/all_sequences/final_edited_pos_mtx.rds")

# Import Motif clust info
mot_clust <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/jeff_non-redundant_pwms.rds")

# Import Motif position
dat <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/sequence_design/Leididi_design/heart/all_sequences/all_batch_motif_pos_ledidi_designed.rds")
setnames(dat, c("seqnames", "motif", "start", "end", "width", "strand", "score"))
dat <- dat[seqnames %in% rownames(edit)]
dat[mot_clust, motif_clust:= i.cluster, on= "motif"]
dat[, seq.idx:= .GRP, seqnames]
dat[, seqnames:= paste0(seqnames, "__", seq.idx)]

# Collapse motifs per cluster
coll <- dat[, collapseBed(.SD), motif_clust, .SDcols= c("seqnames", "start", "end")]
coll[, seq.idx:= tstrsplit(seqnames, "__", keep = 2, type.convert = T)]
coll[, count:= .N, motif_clust]

# Look at motif counts
pl <- unique(coll[, .(count, motif_clust)])[order(count)]
vl_barplot(pl$count, names.arg = pl$motif_clust, ylab= "Create motif counts")
abline(h= 15)

# Select frequent motifs
coll <- coll[count>=15]
setorderv(coll, "count", -1)
coll[, col:= rainbow(.NGRP)[.GRP], motif_clust]
setorderv(coll, "seq.idx")

# Prepare plotting
pl <- coll[, {
  .c <- .SD[, {
    .c <- seq(1, 1001)
    .c <- ifelse(.c %in% unlist(.SD[, start:end, .(start, end)]), col, NA)
    .c[edit[seq.idx,]==0] <- NA
    .c
  }, .(seqnames, seq.idx)]
  .(mat= .(do.call(rbind, split(.c$V1, .c$seqnames))))
}, .(motif_clust, col)]

nSeq <- length(unique(coll$seqnames))
pdf("pdf/motifs_per_sequence.pdf", 40)
vl_par(mgp= c(2, .5, 0),
       mai= c(.9, 2.5, .9, 2))
vl_plot(NA,
        type= "n",
        xlim= c(0, 1001),
        ylim= c(0, nSeq),
        xaxs= "i",
        yaxs= "i",
        frame= F,
        yaxt= "n",
        xlab= "Position",
        ylab= NA)
axis(2, seq(0.5, nSeq-0.5, 1), labels = rev(unique(coll$seqnames)))
pl[, {
  rasterImage(mat[[1]],
              xleft= 0,
              xright = 1001,
              ybottom = 0,
              ytop = nSeq,
              interpolate = F)
}, motif_clust]
segments(0, 1:nSeq, 1001, 1:nSeq)
vl_legend(x= par("usr")[2],
          fill= pl$col,
          legend= pl$motif_clust)
dev.off()