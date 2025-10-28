setwd("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import prediction scores ----
seq.info <- readRDS("Rdata/subbrain_ledidi_design/all_merged_seq_info_old_new.rds")

# Add SOX3 motif counts ----
mot <- readRDS("Rdata/subbrain_ledidi_design/all_merged_motifs_counts_old_new.rds")
seq.info[, SOX3:= mot$SOX3_MOUSE.H11MO.0.C]

# Best design and tissues of interest ----
dat <- seq.info[label %in% c("Design_26_midbrain", "Design_26_forebrain")]
dat[, ID:= .GRP, id]
dat <- melt(
  dat,
  id.vars = c("label", "ID", "SOX3"),
  measure.vars = list(c("act_forebrain", "act_hindbrain", "act_midbrain", "act_neuralTube"),
                      c("acc_forebrain", "acc_hindbrain", "acc_midbrain", "acc_neuralTube")),
  value.name = c("act", "acc")
)
dat[, target_tissue:= tstrsplit(label, "_", keep= 3)]
dat[, model_tissue:= c("forebrain", "hindbrain", "midbrain", "neuralTube")[variable]]

# Retrieve max predicted off-target value ----
dat[, max.other.tissue:= max(act[target_tissue!=model_tissue]), .(label, ID)]

# Compute minimum delta (predicted activity in target tissue - maximum off-target predicted value) ----
dat[, delta:= act-max.other.tissue]

vl_par(mfrow= c(2, 3))
vl_boxplot(act~SOX3,
           dat[target_tissue=="midbrain" & model_tissue=="midbrain"],
           violin= T,
           outline= T,
           xlab= "SOX3 counts",
           ylab= "Predicted activity",
           main= "Midbrain design activity",
           ylim= c(-3, 10))
vl_boxplot(delta~SOX3,
           dat[target_tissue=="midbrain" & model_tissue=="midbrain"],
           violin= T,
           outline= T,
           xlab= "SOX3 counts",
           ylab= "Minimum delta",
           main= "Midbrain design delta\n(vs. strongest off-target)",
           ylim= c(-3, 10))
vl_boxplot(max.other.tissue~SOX3,
           dat[target_tissue=="midbrain" & model_tissue=="midbrain"],
           violin= T,
           outline= T,
           xlab= "SOX3 counts",
           ylab= "Maximum off-target pred. value",
           main= "Midbrain design max off-target",
           ylim= c(-3, 10))
vl_boxplot(act~SOX3,
           dat[target_tissue=="forebrain" & model_tissue=="forebrain"],
           violin= T,
           outline= T,
           xlab= "SOX3 counts",
           ylab= "Predicted activity",
           main= "Forebrain",
           ylim= c(-3, 10))
vl_boxplot(delta~SOX3,
           dat[target_tissue=="forebrain" & model_tissue=="forebrain"],
           violin= T,
           outline= T,
           xlab= "SOX3 counts",
           ylab= "Minimum delta",
           main= "Forebrain design delta\n(vs. strongest off-target)",
           ylim= c(-3, 10))
vl_boxplot(max.other.tissue~SOX3,
           dat[target_tissue=="forebrain" & model_tissue=="forebrain"],
           violin= T,
           outline= T,
           xlab= "SOX3 counts",
           ylab= "Maximum off-target pred. value",
           main= "Forebrain design max off-target",
           ylim= c(-3, 10))

dat[target_tissue=="midbrain" & model_tissue=="midbrain" & SOX3==0][order(delta, decreasing = T)][, .(label, act, max.off.target= max.other.tissue, delta)]
dat[target_tissue=="forebrain" & model_tissue=="forebrain" & SOX3==0][order(delta, decreasing = T)][, .(label, act, max.off.target= max.other.tissue, delta)]
