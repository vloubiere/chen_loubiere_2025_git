setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)
require(GenomicRanges)
require(rtracklayer)

# Import metadata ----
meta <- readxl::read_xlsx("Rdata/metadata_ATACSeq.xlsx")
meta <- as.data.table(meta)[dataset=="bulkENCODE"]

# Import peaks ----
dat <- melt(meta, "tissue", patterns("peaks"))
dat <- dat[, as.data.table(rtracklayer::import(value)), .(tissue, variable)]

# Select canonical chromosomes ----
dat <- dat[seqnames %in% paste0("chr", c(1:19, "X"))]

# Only retain peaks found in both reps and merge ----
conf <- dat[, {
  conf <- vl_intersectBed(.SD[variable=="peaks_merge"],
                          .SD[variable=="peaks_rep1"])
  conf <- vl_intersectBed(conf,
                          .SD[variable=="peaks_rep2"])
  conf
}, tissue]

# Collapse peaks ----
conf[, idx:= vl_collapseBed(.SD, return.idx.only = T)]
conf[, mean_summit:= round(mean(peak)), idx]
conf[, c("start", "end"):= .(min(start), max(end)), idx]
coll <- dcast(conf,
              seqnames+start+end+mean_summit~tissue,
              value.var = "idx",
              fun.aggregate = any)

# Resize (1001bp centered on summit) ----
coll[, start:= start+mean_summit]
res <- vl_resizeBed(coll, "start", upstream = 500, downstream = 500, genome = "mm10")

# Compute coverage and remove strong outliers ----
cov <- meta[, cbind(res[, .(seqnames, start, end)],
                    data.table(score= vl_bw_coverage(res, bw))), tissue]
cov[, norm:= log2(score/sum(score)*1e6), tissue]
cov[, cutoff:= ceiling(quantile(norm, .9995)), tissue]
cov[, check:= all(norm<=cutoff), .(seqnames, start, end)]
cov <- dcast(cov,
             seqnames+start+end+check~tissue,
             value.var = "norm")
print(paste(sum(!cov$check), "strong outlier were remove /", formatC(nrow(cov), big.mark = ","), "total"))
cov <- cov[(check), !"check"]

# Merge peaks and cov ----
dat <- merge(res,
             cov,
             by= c("seqnames", "start", "end"),
             suffixes= c("", ".cov"))

# Clustering ----
layers <- list(peaks= apply(as.matrix(dat[, forebrain:neuralTube]), 2, as.integer),
               cov= as.matrix(dat[, forebrain.cov:neuralTube.cov]))
grid <- somgrid(4, 
                5, 
                "hexagonal", 
                toroidal= T)
init <- lapply(layers, function(x)
{
  set.seed(1)
  x <- x[sample(nrow(x), grid$xdim*grid$ydim), , drop= F]
  return(x)
})
som <- supersom(data = layers, 
                grid= grid,
                init = init,
                user.weights= c(1,5),
                maxNA.fraction = 1)
# Plot
pdf("pdf/ATAC_peaks_clustering_raw.pdf",
    width= 6)
vl_par(mai= c(.9, 2.5,.9,2))
vl_heatmap(dat[, forebrain.cov:neuralTube.cov], 
           row.clusters= som$unit.classif,
           cluster.cols = F,
           show.rownames = F,
           breaks = c(0, 6),
           col= c("blue", "yellow"),
           row.clusters.pos = "left",
           tilt.colnames = T,
           legend.title = "ATAC-Seq lvl (log2)",
           main= "Raw clusters")
dev.off()

# Define classes ----
class <- c("globallyOpen",
           NA,
           "globallyOpen",
           "globallyOpen",
           "globallyOpen",
           "neuralTube",
           "panNeuronal",
           "panNeuronal",
           "heart",
           "hindbrain",
           "hindbrainMidbrainNeuralTube",
           "globallyOpen",
           "limb",
           "forebrain",
           "midbrain",
           NA,
           "heart",
           "forebrain",
           "forebrainHindbrainMidbrain",
           "globallyOpen")
dat$class <- class[som$unit.classif]
dat[, class:= factor(class,
                     c(
                       "heart",
                       "limb",
                       "forebrain",
                       "midbrain",
                       "hindbrain",
                       "neuralTube",
                       "hindbrainMidbrainNeuralTube",
                       "forebrainHindbrainMidbrain",
                       "panNeuronal",
                       "globallyOpen"
                     )
)]
dat[, name:= paste0(class, " (n=", formatC(.N, big.mark = ","), ")"), class]
dat[, name:= factor(name,
                    unique(name[order(class)]))]

# Plot clustering ----
pdf("pdf/ATAC_peaks_clustering.pdf",
    width= 14)
vl_par(mai= c(.9, 2.5,.9,2),
       mfrow= c(1, 2))
vl_heatmap(dat[, forebrain.cov:neuralTube.cov], 
           row.clusters= paste0(dat$class, ".", som$unit.classif),
           cluster.cols = F,
           show.rownames = F,
           breaks = c(0, 6),
           col= c("blue", "yellow"),
           row.clusters.pos = "left",
           tilt.colnames = T,
           legend.title = "ATAC-Seq lvl (log2)",
           main= "Full clustering")

# Plot refined clusters ----
vl_heatmap(dat[, .(heart.cov, limb.cov, forebrain.cov, midbrain.cov, hindbrain.cov, neuralTube.cov)],
           row.clusters= dat$name,
           cluster.cols = F,
           show.rownames = F,
           breaks = c(0, 6),
           col= c("blue", "yellow"),
           row.clusters.pos = "left",
           tilt.colnames = T,
           legend.title = "ATAC-Seq lvl (log2)",
           main= "Simplified clusters")
dev.off()
print(paste(formatC(sum(is.na(dat$class)), big.mark = ","), "regions / ", formatC(nrow(dat), big.mark = ","), "not classified -> removed"))

# Update overlaps ----
dat[, hindbrain:= as.integer(class %in% c("globallyOpen", "panNeuronal") | grepl("hindbrain", class, ignore.case = T))]
dat[, neuralTube:= as.integer(class %in% c("globallyOpen", "panNeuronal") | grepl("neuralTube", class, ignore.case = T))]
dat[, midbrain:= as.integer(class %in% c("globallyOpen", "panNeuronal") | grepl("midbrain", class, ignore.case = T))]
dat[, forebrain:= as.integer(class %in% c("globallyOpen", "panNeuronal") | grepl("forebrain", class, ignore.case = T))]
dat[, heart:= as.integer(class=="globallyOpen" | grepl("heart", class, ignore.case = T))]
dat[, limb:= as.integer(class=="globallyOpen" | grepl("limb", class, ignore.case = T))]

# Save ----
setorderv(dat, c("seqnames", "start", "end"))
saveRDS(dat[!is.na(class), .(peakID= paste0(seqnames, ":", start, "-", end), seqnames, start, end, class, heart, limb, forebrain, midbrain, hindbrain, neuralTube)],
        "db/peaks/bulkENCODE_confident_ATAC_peaks.rds")
