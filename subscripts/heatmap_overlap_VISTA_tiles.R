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
# cols <- c("heart", "limb", "FacialMes", "Nose", "DRG", "forebrain", "midbrain", "hindbrain", "neuralTube")
ov <- dat[, {
  m <- as.matrix(.SD)+0
  t(m)%*%m
}, .SDcols= cols]

# Overlap ----
vl_heatmap(
  ov,
  breaks= seq(0, 400, length.out= 21),
  cluster.rows = F,
  show.numbers = T,
  numbers.cex = 0.7*5/8,
  legend.title = "Active VISTA tile",
  pdf.file = "pdf/0_paper/heatmap_overlap_VISTA_tiles.pdf",
  pdf.cell.size = .2
)
