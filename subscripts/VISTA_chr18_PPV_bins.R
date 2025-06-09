setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Randomly sample chromosome 18
rdm <- GenomeInfoDb::seqlengths(BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10)
rdm <- as.data.table(rdm, keep.rownames= "seqnames")
rdm <- rdm[seqnames=="chr18"]
set.seed(1)
rdm <- rdm[, .(start= sample(rdm-1001, 1e6)), seqnames]
rdm[, end:= start+1000]
set.seed(1)
rdm[, strand:= sample(c("+", "-"), 1e6, replace = T)]

# Remove centromeres and telomeres  ----
heterochromatin <- fread("/groups/stark/shenzhi.chen/db/blacklisted/telomere_centromeres_mm10.bed",
                         col.names =  c("seqnames", "start", "end", "name"))
heterochromatin[, start:= start-10e3]
heterochromatin[, end:= end+10e3]
rdm <- vl_intersectBed(rdm,
                       heterochromatin,
                       invert = T)

# Remove overlaps with VISTA active tiles ----
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
vista[, start:= start-100]# Extended to reflect later augmentation
vista[, end:= end+100]# Extended to reflect later augmentation
rdm <- vl_intersectBed(rdm,
                       vista[!is.na(start)],
                       invert = T)

# Remove open peaks
ATAC <- readRDS("db/peaks/bulkENCODE_confident_ATAC_peaks.rds")
ATAC <- vl_resizeBed(ATAC, "center", 1400, 1400, genome = "mm10") # Extended to reflect later augmentation
rdm <- vl_intersectBed(rdm,
                       ATAC,
                       invert = T)

# Randomly sample 300000
set.seed(1)
rdm <- rdm[sample(.N, 3e5)]

# Extract sequences 
rdm[, seq:= vl_getSequence(rdm, "mm10")]
rdm[, ID:= paste0(seqnames, ":", start, "-", end, ":", strand, "__inactive")]
rdm[, score:= 0]

# Save fasta file
seqinr::write.fasta(sequences = as.list(rdm[, seq]), 
                    names = rdm[, ID], 
                    file.out = "db/fasta/bulkATAC/VISTA/PPV_bins_chr18.fa")

# Save fasta file
rdm[, ID:= paste0(seqnames, ":", start, "-", end, ":", strand, "__inactive")]
seqinr::write.fasta(sequences = as.list(rdm[, seq]), 
                    names = rdm[, ID], 
                    file.out = "db/fasta/bulkATAC/VISTA/PPV_bins_chr18.fa")

# Save observed file
fwrite(rdm[, .(class= "inactive", ID, score)],
       "db/observed/bulkATAC/VISTA/PPV_bins_chr18.txt",
       sep= "\t",
       na = NA)

# Save bed file
rtracklayer::export(GRanges(rdm[, .(seqnames, start, end, strand, name= ID, score)]),
                    "db/bed/bulkATAC/VISTA/PPV_bins_chr18.bed")  
