setwd("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import prediction scores ----
seq.info <- readRDS("Rdata/subbrain_ledidi_design/all_sequence_information.rds")

# Import prediction scores ----
dat <- melt(
  seq.info,
  id.vars = c(
    "setting", # LEDIDI SETTINGS (see excel file)
    "target_tissue", # Target active tissue
    "model_tissue", # Model used to predict the activity of the sequence
    "fold",
    "rep",
    "location", # unique sequence ID (string containing all LEDIDI edits)
    "sequences" # DNA sequence
    ),
  measure.vars = c("accessibility", "sigmoid_activity"),
)

# Mean prediction per Design/tissue ----
dat <- dat[, .(mean_prediction = mean(value, na.rm = TRUE)),
           by = .(variable, setting, target_tissue, model_tissue)]

# Retrieve max predicted off-target value ----
dat[, max.other.tissue:= max(mean_prediction[model_tissue!=target_tissue]), .(variable, setting, target_tissue)]

# Compute minimum delta (predicted activity in target tissue - maximum off-target predicted value) ----
dat[, delta:= mean_prediction-max.other.tissue]

# Dcast matrices ----
# Activity
act.delta <- dcast(dat[target_tissue==model_tissue & variable=="sigmoid_activity"], setting~target_tissue, value.var = "delta")
act.delta <- as.matrix(act.delta, 1)
pred.act <- dcast(dat[target_tissue==model_tissue & variable=="sigmoid_activity"], setting~target_tissue, value.var = "mean_prediction")
pred.act <- as.matrix(pred.act, 1)
# Accessibility
acc.delta <- dcast(dat[target_tissue==model_tissue & variable=="accessibility"], setting~target_tissue, value.var = "delta")
acc.delta <- as.matrix(acc.delta, 1)
pred.acc <- dcast(dat[target_tissue==model_tissue & variable=="accessibility"], setting~target_tissue, value.var = "mean_prediction")
pred.acc <- as.matrix(pred.acc, 1)

# Plot heatmaps ----
pdf("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/pdf/0_paper/heatmap_mean_minimum_delta_sub_CNS_design.pdf", width = 9)
# Activity
vl_par(mfrow= c(1,3), omi= c(0,0,0,1))
vl_heatmap(act.delta, cluster.rows = scale(act.delta), show.numbers = round(act.delta, 1), legend.title = "Min. act. delta")
vl_heatmap(scale(act.delta), show.numbers = round(act.delta, 1), legend.title = "Scaled\nMin. act. delta")
vl_heatmap(pred.act, cluster.rows = scale(act.delta), show.numbers = round(pred.act, 1), legend.title = "Mean pred. act.")
# Accessibility
vl_heatmap(acc.delta, cluster.rows = scale(acc.delta), show.numbers = round(acc.delta, 1), legend.title = "Min. acc. delta")
vl_heatmap(scale(acc.delta), show.numbers = round(acc.delta, 1), legend.title = "Scaled\nMin. acc. delta")
vl_heatmap(pred.acc, cluster.rows = scale(acc.delta), show.numbers = round(pred.acc, 1), legend.title = "Mean pred. acc.")
dev.off()