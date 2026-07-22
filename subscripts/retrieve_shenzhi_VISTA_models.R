setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)


# "/groups/stark/shenzhi.chen/model/deepstarr2-mouse/Data/Input/single_task/Muti_tissue_accessibility_model_20240617/Heart/Enhancer_model_shared_20240805/Results_fold01_score_rep1/fold01_sequences_test.fa_predictions_Model.txt"

# Retrieve test set predictions ATAC-Seq models ----
meta <- data.table(tissue= c("heart", "limb", "neuralTube", "midbrain", "hindbrain", "forebrain"),
                   Tissue= c("Heart", "Limb", "Neuraltube", "Midbrain", "Hindbrain", "Forebrain"))
meta[, folder:= paste0("/groups/stark/shenzhi.chen/model/deepstarr2-mouse/Data/Input/single_task/Muti_tissue_accessibility_model_20240617/", Tissue, "/Enhancer_model_shared_20240805")]
meta <- meta[, .(predFile= list.files(folder, "_sequences_test.fa_predictions_Model.txt$", full.names = T, recursive= T)), (meta)]

# Retrieve conditions ----
meta[, info:= tstrsplit(predFile, "/", keep= 13)]
meta <- meta[!grepl("heart$", info)]
meta[, c("fold", "rep"):= tstrsplit(info, "_", keep= c(2,4))]
meta[grepl("aug_1$", info), aug:= "aug1"]
meta[grepl("aug_2$", info), aug:= "aug2"]
meta[grepl("aug_3$", info), aug:= "aug3"]
meta[is.na(aug), aug:= "No"]

# Retrieve test set observed values ----
meta[, obsFile:= paste0("db/scores/VISTA/", tissue, "/", fold, "_sequences_activity_test.txt"), .(tissue, fold)]

# Import observed and predicted and merge ----
obs <- meta[, fread(obsFile), .(tissue, fold, obsFile)]
pred <- meta[, fread(predFile), .(tissue, fold, rep, aug, predFile)]
dat <- merge(obs[, .(tissue, fold, class, ID= nameID, active= label=="Active")],
             pred[, .(tissue, fold, rep, aug, ID= location, predicted= Predictions)],
             by= c("tissue", "fold", "ID"),
             allow.cartesian= T)

# Save data ----
saveRDS(dat,
        "db/predictions/20240807_VISTA_model_testSet_predictions.rds")

# Retrieve shared test set predictions ----
shared <- data.table::copy(meta)
shared[, predFile:= gsub("_sequences_test.fa_predictions_Model.txt$", "_sequences_sharedTest.fa_predictions_Model.txt", predFile)]
shared[, obsFile:= gsub("_sequences_activity_test.txt$", "_sequences_activity_sharedTest.txt", obsFile)]

# Import observed and predicted and merge ----
obs <- shared[, fread(obsFile), .(tissue, fold, obsFile)]
pred <- shared[, fread(predFile), .(tissue, fold, rep, aug, predFile)]
dat <- merge(obs[, .(tissue, fold, class, ID= nameID, active= label=="Active")],
             pred[, .(tissue, fold, rep, aug, ID= location, predicted= Predictions)],
             by= c("tissue", "fold", "ID"),
             allow.cartesian= T)

# Save data ----
saveRDS(dat,
        "db/predictions/20240807_VISTA_model_sharedTest_predictions.rds")