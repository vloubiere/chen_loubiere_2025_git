setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Retrieve test set predictions ATAC-Seq models ----
meta <- data.table(tissue= c("heart", "limb", "neuralTube", "midbrain", "hindbrain", "forebrain"),
                   Tissue= c("Heart", "Limb", "Neuraltube", "Midbrain", "Hindbrain", "Forebrain"))
meta[, folder:= paste0("/groups/stark/shenzhi.chen/model/deepstarr2-mouse/Data/Input/single_task/Muti_tissue_accessibility_model_20240617/", Tissue, "/data_output_merged_2_weighted_20240802")]
meta <- meta[, .(predFile= list.files(folder, "_sequences_test.fa_predictions_enhancer_Model.txt$", full.names = T, recursive= T)), (meta)]

# Retrieve conditions ----
meta[, fold:= tstrsplit(predFile, "/", keep= 13)]
meta[, c("fold", "rep"):= tstrsplit(fold, "_", keep= c(2,5))]

# Retrieve test set observed values ----
meta[, obsFile:= paste0("db/scores/ATAC/", tissue, "/", fold, "_sequences_activity_test.txt"), .(tissue, fold)]

# Import observed and predicted and merge ----
obs <- meta[, fread(obsFile), .(tissue, fold, obsFile)]
pred <- meta[, fread(predFile), .(tissue, fold, rep, predFile)]
dat <- merge(obs[, .(tissue, fold, class, ID, observed= score)],
             pred[, .(tissue, fold, rep, ID= location, predicted= Predictions)],
             by= c("tissue", "fold", "ID"),
             allow.cartesian= T)

# Save data ----
saveRDS(dat,
        "db/predictions/20240807_ATACSeq_model_testSet_predictions.rds")

# Retrieve test chromosome predictions ----
chr <- meta[fold=="fold07"]
chr[Tissue=="Neuraltube", Tissue:= "Neural-tube"]
chr[, obsFile:= paste0("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/observed/Chr11_", Tissue, "_observed_act.txt"), Tissue]
chr[, predChr11:= list.files(dirname(predFile), "^whole_chr11", full.names = T), predFile]

# Import observed predicted and save ----
obs <- chr[, fread(obsFile), .(tissue, fold, obsFile)]
pred <- chr[, fread(predFile), .(tissue, fold, rep, predFile= predChr11)]
dat2 <- merge(obs[, .(tissue, fold, ID, observed= score)],
              pred[, .(tissue, fold, rep, ID= location, predicted= Predictions)],
              by= c("tissue", "fold", "ID"),
              allow.cartesian= T)

# Save data ----
saveRDS(dat2,
        "db/predictions/20240807_ATACSeq_model_chr11_predictions.rds")