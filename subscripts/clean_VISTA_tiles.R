setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
source("git_deepATAC/function/augmentation_function_tiling_sliding_window.R")
require(vlfunctions)

# Import vista tiles ----
vista <- readxl::read_excel("/groups/stark/shenzhi.chen/db/VISTA_enhancer_dataset/VISTA2024_AllTissuesReferenceAlleles.xlsx")
vista <- as.data.table(vista)
vista <- vista[, .(class= "vistaTile",
                   peakID= Vista.ID,
                   genome= fcase(grepl("^hs", Vista.ID), "hg38",
                                 grepl("^mm", Vista.ID), "mm10"),
                   coor_hg38= Coordinates_hg38,
                   coor_mm10= Coordinates.mm10,
                   heart= Heart,
                   limb= Limb,
                   forebrain= Forebrain,
                   midbrain= Midbrain,
                   hindbrain= Hindbrain,
                   neuralTube= NeuralTube)]
vista <- vista[!is.na(genome)]

# Remove mutated enhancers ----
mut <- fread("db/public/VISTA_experiments_20250609.tsv")[grepl("mutagenesis", description), vista_id]
vista <- vista[!(peakID %in% mut)]

# Extract genomic coordinates ----
vista[genome=="mm10", c("seqnames", "start", "end"):= importBed(coor_mm10)[, 1:3]]
vista[genome=="hg38", c("seqnames", "start", "end"):= importBed(coor_hg38)[, 1:3]]

# Remove long enhancers (>5kb) ----
vista[, width:= end-start+1]
vista <- vista[width<=5000]

# Resize short enhancers to 1.5kb ----
# We will be using 1001bp tiles for augmentation, so center +/-250 will always be there
vista[width<1500, start:= round((start+end)/2)-750]
vista[width<1500, end:= start+1501]
vista[, width:= NULL]

# Order and save ----
setcolorder(vista,
            c("class", "genome", "peakID", "seqnames", "start", 'end'))
saveRDS(vista,
        "db/peaks/vista_tiles_clean.rds")



