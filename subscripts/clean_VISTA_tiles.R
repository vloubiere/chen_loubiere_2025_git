setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
source("git_deepATAC/function/augmentation_function_tiling_sliding_window.R")
require(vlfunctions)

# Import vista tiles ----
vista <- readxl::read_excel("/groups/stark/shenzhi.chen/db/VISTA_enhancer_dataset/VISTA2024_AllTissuesReferenceAlleles.xlsx")
vista <- as.data.table(vista)
vista[, class:= "vistaTile"]

# Remove mutated enhancers ----
mut <- fread("/groups/stark/shenzhi.chen/db/VISTA_enhancer_dataset/download.csv")
vista <- vista[!(Vista.ID %in% mut[grepl("mutagenesis", `Element Description`), `Vista ID`])]

# Extract mm10 genomic coordinates and human ones ----
vista[, genome:= fcase(grepl("^hs", Vista.ID), "hg38",
                       grepl("^mm", Vista.ID), "mm10")]
vista[, c("seqnames", "start", "end"):= tstrsplit(Coordinates.mm10, ":|-", type.convert = T)]
vista[, c("hs.seqnames", "hs.start", "hs.end"):= tstrsplit(Coordinates_hg38, ":|-", type.convert = T)]

# Remove long enhancers (>5kb) ----
vista[genome=="mm10", width:= end-start+1]
vista[genome=="hg38", width:= hs.end-hs.start+1]
vista <- vista[width<5000]

# Resize short enhancers to 1.5kb ----
# With 1001bp tiles, center +/-250will be there in all augmented sequences
mm.short <- vista[, end-start+1]
vista[mm.short<1500, start:= round((start+end)/2)-750]
vista[mm.short<1500, end:= start+1501]
hs.short <- vista[, hs.end-hs.start+1]
vista[hs.short<1500, hs.start:= round((hs.start+hs.end)/2)-750]
vista[hs.short<1500, hs.end:= hs.start+1501]

# Final coordinates ----
vista[genome=="mm10", coor:= paste0(seqnames, ":", start, "-", end)]
vista[genome=="hg38", coor:= paste0(hs.seqnames, ":", hs.start, "-", hs.end)]

# Save ----
saveRDS(vista[, .(peakID= Vista.ID,
                  genome, seqnames, start, end, class, coor,
                  heart= Heart, limb= Limb, forebrain= Forebrain, midbrain= Midbrain, hindbrain= Hindbrain, neuralTube= NeuralTube)],
        "db/peaks/vista_tiles_clean.rds")



