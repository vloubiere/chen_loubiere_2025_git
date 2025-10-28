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

# Upset plot
pdf("pdf/0_paper/upset_plot_midbrain.pdf", 4.5, 5)
vl_par(mai= c(3,2,.9,.9))
.m <- melt(dat[midbrain==1, .(peakID, forebrain, hindbrain, neuralTube, midbrain)],
             id.vars = "peakID")
.m <- .m[value==1]
upsetPlot(split(.m$peakID, .m$variable), cex.grid = .7, grid.hex = .7)
dev.off()

# Pie chart
midbrain <- dat[midbrain==1]
pdf("pdf/0_paper/pie_chart_midbrain_enhancers.pdf", 4.5, 5)
vl_par(mai= rep(1.5, 4))
specific <- fcase(
  midbrain[, midbrain+forebrain+hindbrain+neuralTube]==1, "Midbrain-specific",
  midbrain[, midbrain+forebrain+hindbrain+neuralTube]==4, "Pan-CNS",
  midbrain[, midbrain+forebrain+hindbrain+neuralTube]>1, "Shared with >= 1 CNS sub-region"
)
specific <- factor(specific, c("Midbrain-specific", "Shared with >= 1 CNS sub-region", "Pan-CNS"))
vl_pie(specific,
       init.angle= 180, 
       labels = "p", 
       clockwise= T,
       col= adjustcolor(c("limegreen", "lightskyblue", "tomato"), .3))
dev.off()