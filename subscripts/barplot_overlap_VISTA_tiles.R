setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import data ----
dat <- readRDS("db/folds/bulkATAC_folds.rds")
dat <- dat[group=="vista"]

# Add missing labels ----
add <- readxl::read_xlsx("/groups/stark/shenzhi.chen/db/VISTA_enhancer_dataset/VISTA2024_AllTissuesReferenceAlleles.xlsx")
add <- as.data.table(add)
dat <- merge(dat, add[, .(peakID= Vista.ID, FacialMes, Nose, DRG)], all.x= T)

# Overlap ----
cols <- c("heart", "limb", "forebrain", "midbrain", "hindbrain", "neuralTube")
ov <- as.matrix(dat[, lapply(.SD, sum), .SDcols= cols])

# Overlap ----
pdf("pdf/_revision/barplot_number_VISTA_enh_per_tissue.pdf", 3.25, 3.25)
vl_par()
vl_barplot(ov, col= "lightgrey", border= NA, width = .8, ylab= "Active VISTA sequences")
dev.off()