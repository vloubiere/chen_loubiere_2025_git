setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[set %in% c("test", "validation", "training")]
meta <- meta[dataset %in% "activity"]
meta <- meta[tissue %in% c("heart", "midbrain", "limb")]
meta <- meta[ID=="model1_bulkATAC_tsx3Aug_2xBal_noW"]

# ChIPSeq peaks ----
peaks <- readRDS("Rdata/revision_ChIPseq_peaks.rds")
peaks <- peaks[signalValue>5]

# ATAC peaks ----
ATAC <- readRDS("db/peaks/bulkENCODE_confident_ATAC_peaks.rds")
ATAC <- melt(ATAC, c("seqnames", "start", "end"), c("limb", "midbrain", "heart"))
ATAC <- ATAC[value==1]

# Promoters ----
if(!exists("ext_prom")) {
  gtf <- rtracklayer::import("../../genomes/Mus_musculus/GENCODE/gencode.vM25.annotation.gtf.gz")
  gtf <- as.data.table(gtf)
  gtf <- gtf[type=="transcript"]
  ext_prom <- resizeBed(gtf[, .(seqnames, start, end)], "start", 10000, 10000)
  ext_prom <- unique(ext_prom)
}

# Tests ----
PPV <- meta[, {
  
  # Import bed file ----
  .c <- rbindlist(lapply(unique(bed_file), importBed))
  .c[, c("ID", "strand.1", "tile", "score"):= tstrsplit(name, "_|__|:")]
  .c[, score:= as.numeric(ifelse(score=="active", 1, 0))]
  
  # Compute overlaps with enhancer marks ----
  .c$K27Ac_ov <- covBed(.c, peaks[.BY, on= "tissue"][cdition=="H3K27ac"])>0
  .c$K4me1_ov <- covBed(.c, peaks[.BY, on= "tissue"][cdition=="H3K4me1"])>0
  .c$K4me3_ov <- covBed(.c, peaks[.BY, on= "tissue"][cdition=="H3K4me3"])>0
  
  # Compute overlaps with genomic annotations ----
  .c$ext_prom <- covBed(.c, ext_prom)>0
  .c$CTCF_ov <- covBed(.c, peaks[cdition=="CTCF"])>0 # Union of all tissues
  
  # Tissue-specific ATAC-seq peaks ----
  .c$ATAC_ov <- covBed(.c, ATAC[variable==tissue])>0
  .c$ATAC_other_tissue <- covBed(.c, ATAC[variable!=tissue])>0
  
  # Compute TL enhancer labels ----
  .c[, TL0_all:= ATAC_ov] # Accessible in the tissue
  .c[, TL1_K27Ac:= TL0_all & K27Ac_ov]
  .c[, TL2_K27Ac_K4me1:= TL0_all & K27Ac_ov & K4me1_ov]
  .c[, TL3_K27Ac_K4me1_noK4me3:= TL0_all & K27Ac_ov & K4me1_ov & !K4me3_ov]
  .c[, TL4_prom10kb:= TL0_all & (!ext_prom)]
  .c[, TL5_prom10kb_noCTCF:= TL0_all & (!ext_prom) & (!CTCF_ov)]
  .c[, TL6_tsATAC:= TL0_all & (!ATAC_other_tissue)]
  .c[, TL7_HTM_tsATAC:= TL3_K27Ac_K4me1_noK4me3 & TL6_tsATAC] # Combine HTMs and tissue-specific ATAC
  
  # Reformat ----
  var <- melt(.c, id.vars= c("ID", "score"), patterns("TL"), value.name = "score.new")
  var[, score:= as.logical(score)]
  
  # Aggregate ----
  agg <- var[, .(score= any(score), score.new= any(score.new)), .(ID, variable)]
  
  # # Number of positive ----
  # agg[, variable:= gsub("TL0_|TL1_|TL2_|TL3_|TL4_|TL5_|TL6_|TL7_", "", variable)]
  # agg[, TL:= paste0("TL", 8:14)[.GRP], variable]
  # agg[, c("TRAINING_POS", "TRAINING_NEG"):= {
  #   .t <- fread(paste0("db/observed/bulkATAC/", TL, "/heart/fold01_sequences_activity_training.txt"))
  #   .(sum(.t$class=="active"), sum(.t$class=="inactive"))
  # }, TL]
  
  # Return ----
  agg[, {
    .(VISTA_POS= sum(score),
      VISTA_NEG= sum(!score),
      intersection_VISTA_all= sum(score.new),
      TRUE_POS= sum(score & score.new),
      FALSE_POS= sum(!score & score.new),
      FALSE_NEG= sum(score & !score.new))
  }, variable]
}, tissue]

# Format
final <- PPV[variable %in% c("TL0_all", "TL5_prom10kb_noCTCF", "TL2_K27Ac_K4me1", "TL6_tsATAC"), .(
  Tissue= tissue,
  VISTA_active= VISTA_POS,
  VISTA_inactive= VISTA_NEG,
  rule= variable,
  TP= TRUE_POS,
  FP= FALSE_POS,
  Recall= TRUE_POS/(TRUE_POS+FALSE_NEG)*100,
  Precision= TRUE_POS/(TRUE_POS+FALSE_POS)*100
)]
final[, rule:= factor(rule, c("TL0_all", "TL5_prom10kb_noCTCF", "TL2_K27Ac_K4me1",  "TL6_tsATAC"))]

# Plot ----
Cc <- c("grey40", "cornflowerblue", "limegreen", "tomato")
pdf("pdf/_revision/barplot_precision_recall_TL.pdf", height = 3)
vl_par(mfrow= c(1,2), las= 1)
barplot(Precision~rule+Tissue, final, col= Cc, ylab= "Precision (%)", beside= T)
barplot(Recall~rule+Tissue, final, col= Cc, ylab= "Recall (%)", beside= T)
vl_legend(legend = levels(final$rule), fill= Cc)
dev.off()

