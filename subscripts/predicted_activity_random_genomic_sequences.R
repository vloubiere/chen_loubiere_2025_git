setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & tissue %in% c("heart", "limb", "midbrain")]
meta <- meta[set=="NegGenomicRegions"]

# Import data ----
dat <- meta[, {
  fread(pred_file)
}, .(pred_file, tissue)]

# Average ----
agg <- dat[, .(pred= mean(Predictions)), .(tissue, location)]
agg[, location:= factor(location, unique(location))]
agg <- dcast(agg, location~tissue, value.var = "pred")
agg[, sequence:= seqinr::read.fasta(meta$fa_file[1], as.string = T)]
agg[, sequence:= toupper(as.character(sequence))]

# Call motifs ----
mot <- readRDS("../../motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")
sel.mot <- c("MEF2A_MOUSE.H11MO.0.A", "TWIST1_MA1123.1", "SOX3_MOUSE.H11MO.0.C")
sel.mot <- mot$pwms_log_odds[match(sel.mot, mot$meta$Motif)]
counts <- vl_motifCounts(
  sequences = agg$sequence,
  pwm_log_odds = sel.mot,
  bg = "genome",
  genome = "mm10",
  p.cutoff = 1e-5
)
counts[, location:= agg$location]

# Count labels ----
res <- melt(counts, id.vars= "location")
res[, class:= {
  label <- as.character(value)
  if(max(value)>2)
    label[value>=2] <- ">=2"
  label <- factor(label, c("0", "1", ">=2"))
  class <- table(label)
  class <- paste0(names(class), " (n= ", formatC(class, big.mark= ","), ")")
  factor(class[label], class)
}, variable]

# Add predicted activity ----
final <- melt(agg, id.vars = "location", measure.vars = c("heart", "limb", "midbrain"), value.name = "act", variable.name = "tissue")
final <- merge(final, res[, .(location, variable, class)], by= "location", allow.cartesian= T)
final <- final[
  (tissue=="heart" & variable=="MEF2__MEF2A_MOUSE.H11MO.0.A") |
    (tissue=="limb" & variable=="Ebox/CAGCTG__TWIST1_MA1123.1") |
    (tissue=="midbrain" & variable=="SOX/1__SOX3_MOUSE.H11MO.0.C")
]

# Plot ----
pdf("pdf/_revision/boxplot_predicted_activity_random_genomic_sequences_motif_counts.pdf", 3.5, 2.75)
vl_par(mfrow= c(1,3), mai= c(.9, .1, .9, .2), omi= c(0, .5, 0, 0), cex.main= 1)
final[, {
  .c <- droplevels(class)
  vl_boxplot(
    act~.c,
    main= paste0(strsplit(as.character(variable), "_")[[1]][3], " in ", tissue[1]),
    tilt.names= T,
    compute.pval= list(c(1,2), c(1,3))
  )
  if(.GRP==1)
    title(ylab= "Activity score", xpd= NA)
  print("")
}, variable]
dev.off()