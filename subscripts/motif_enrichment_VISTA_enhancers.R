setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)
require(dtBedTools)

# Import VISTA sequences ----
vista <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/vista_tiles_clean.rds")
vista$sequences <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/vista_tiles_with_sequences.rds")$sequences
tissues <- c("heart", "limb", "forebrain", "midbrain", "hindbrain", "neuralTube")
vista[, Ntissues:= rowSums(.SD), .SDcols= tissues]
vista <- vista[Ntissues <= 1]
class <- melt(vista, "peakID", tissues)[value==1]
vista[class, class:= i.variable, on= "peakID"]
vista[class=="vistaTile", class:= "neg"]

# Generate random controls ----
rdm <- vl_control_regions_BSgenome(bed = vista[!is.na(seqnames)],
                                   genome = "mm10",
                                   no.overlap = TRUE)
rdm[, sequences:= vl_getSequence(rdm, "mm10")]

# Count motifs ----
motifs <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/jeff_non-redundant_pwms_jaspar_all.rds")
vCount <- "db/counts/motifs/vistaCounts.rds"
if(!file.exists(vCount))
{
  .c <- vl_motif_counts(sequences = vista$sequences,
                        pwm_log_odds = motifs,
                        genome = "mm10",
                        p.cutoff = 1e-4)
  saveRDS(.c,
          vCount)
}
vCount <- readRDS(vCount)
cCount <- "db/counts/motifs/vistaControlCounts.rds"
if(!file.exists(cCount))
{
  .c <- vl_motif_counts(sequences = rdm$sequences,
                        pwm_log_odds = motifs,
                        genome = "mm10",
                        p.cutoff = 1e-4)
  saveRDS(.c,
          cCount)
}
cCount <- readRDS(cCount)

# Compute enrichment ----
enr <- vl_motif_cl_enrich(counts.list = c(split(vCount, vista$class), list(control= cCount)),
                          names = colnames(vCount),
                          control.cl = "control")

# Select motifs with the best padj per cluster ----
cl <- readxl::read_xlsx("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/motif_annotations.xlsx", sheet = 2)
cl <- as.data.table(cl)
enr[cl, cluster:= i.Cluster_ID, on= "variable==Motif"]
enr[, best:= padj==min(padj), cluster]
enr <- enr[variable %in% enr[(best), variable]]

# Plot
mat <- dcast(enr, variable~cl, value.var = "padj")
mat <- as.matrix(mat, 1)

pdf("pdf/motif_enrichment_vista.pdf", height = 30)
vl_par(mai= c(.9, 2, .9, 2), cex.axis = .5)
vl_heatmap(-log10(mat),
           breaks = c(-1, 0, 15),
           display.numbers = TRUE,
           cluster.cols = FALSE, display.numbers.cex = .5)
dev.off()


# Import Designed sequences ----
designed <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/enhancer_design/evolutional_design/6000_designed_enhancers.rds")
designed <- data.table(peakID= tstrsplit(names(designed), "/", keep= 5),
                       seq= as.character(designed))
