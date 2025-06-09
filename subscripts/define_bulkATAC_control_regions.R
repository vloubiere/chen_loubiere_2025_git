setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Bin the whole genome  ----
bins <- vl_binBSgenome(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10,
                       bins.width = 1001,
                       steps.width = 1000)
bins <- vl_resizeBed(bins, "center", 1400, 1400, genome= "mm10") # Account for later augmentation

# Select canonical chromosomes ----
bins <- bins[seqnames %in% paste0("chr", c(1:19, "X"))]

# Remove centromeres and telomeres  ----
heterochromatin <- fread("/groups/stark/shenzhi.chen/db/blacklisted/telomere_centromeres_mm10.bed",
                         col.names =  c("seqnames", "start", "end", "name"))
heterochromatin[, start:= start-10e3]
heterochromatin[, end:= end+10e3]
bins <- vl_intersectBed(bins, heterochromatin, invert = T)

# Remove the bins that are less than 1.5kb bp away from closest peak summit to account for augmentation ----
peaks <- readxl::read_xlsx("Rdata/metadata_ATACSeq.xlsx")
peaks <- as.data.table(peaks)[dataset=="bulkENCODE"]
peaks <- na.omit(unique(unlist(peaks[, peaks_merge:peaks_rep2])))
peaks <- lapply(peaks, function(x) as.data.table(rtracklayer::import(x)))
peaks <- rbindlist(peaks)
peaks[, start:= start+peak] # Summit
peaks <- vl_resizeBed(peaks, "start", 1500, 1500) # Account for later augmentation
bins <- vl_intersectBed(bins,
                        peaks,
                        invert = T)

# Remove bins overlapping VISTA tiles ----
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
bins <- vl_intersectBed(bins,
                        vista[!is.na(seqnames)],
                        invert = T)

# Randomly select 1M negative regions ----
set.seed(1)
bins <- bins[sample(.N, 1e6)]

# Save ----
bins <- vl_resizeBed(bins, "center", 500, 500)
bins[, class:= "globallyClosed"]
setorderv(bins,
          c("seqnames", "start", "end"))
saveRDS(bins[, .(peakID= paste0(seqnames, ":", start, "-", end), seqnames, start, end, class,
                 heart= 0L, limb= 0L, forebrain= 0L, midbrain= 0L, hindbrain= 0L, neuralTube= 0L)],
        "db/peaks/bulkENCODE_control_regions.rds")
