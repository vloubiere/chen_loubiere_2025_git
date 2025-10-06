# Set workpath ----
setwd("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/clean_version/")
library(Biostrings)
devtools::load_all("/groups/stark/vloubiere/vlite/")
options(datatable.prettyprint.char = 20)

# Import metadata ----
meta <- readRDS("Rdata/data_preprocess/3.metadata_predicted_accessibility.rds")

# In the end we only used 3 folds and 2 replicates ----
meta <- meta[fold %in% c("fold01","fold02","fold03") & replicate %in% c("rep1", "rep2")]

# Remove useless columns ----
meta$dataset <- meta$Batch <-meta$augmentation <- meta$augmentation.detail <- meta$balancing <- meta$balancing.detail <- meta$weight <- meta$weight.detail <- NULL

# Add missing observed values for the full chr18 (used as test set) ---
meta[, obs_access_chr18 := paste0("db/testing_dataset/chr18/", tissue, "_chr18_bins_access_obs.txt")]
meta[, fa_access_chr18 := paste0("db/testing_dataset/chr18/", tissue, "_chr18_bins_access_obs.fa")]
meta[, bed_access_chr18 := paste0("db/testing_dataset/chr18/", tissue, "_chr18_bins_access_obs.bed")]

# Add fasta files containing random sequences (used as extra test set) ----
meta[, fa_access_random:= "db/testing_dataset/random_sequences_600k.fasta"]
meta[, fa_VISTA_random:= "db/testing_dataset/random_sequences_600k.fasta"]

# Add missing VISTA predictions ----
meta[, c("pred_VISTA_test", "pred_VISTA_random"):= {
  # folder 
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  .(
    file.path(folder, paste0(fold, "_sequences_test.fa_predictions_enhancer_Model.txt")),
    file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
  )
}, .(ID, tissue, fold, replicate)]

# Add contribution files ----
meta[, contrib_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  file.path(folder, paste0("/Model_", fold,"_sequences_test.fa_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5"))
}]
meta[, contrib_access_test := {
  folder <- paste0("result/model/ATAC_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  file.path(folder, paste0("/Model_", fold ,"_sequences_test.fa_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5"))
}]

# Add activity prediction from accessibility model ----
meta[, predActFromAcc_VISTA_test := {
  folder <- paste0("result/model/ATAC_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "/act")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[, predActFromAcc_VISTA_random := {
  folder <- paste0("result/model/ATAC_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]

# Add activity prediction from randomly initialized models  ----
meta[, predActFromRandomIni_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_init_random")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[, predActFromRandomIni_VISTA_random := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_init_random")
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]

# Add activity prediction from mismatched tissue TL models (MTL) ----
meta[tissue== "heart", predActFromMTL_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_neuralTube")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[tissue== "limb", predActFromMTL_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_neuralTube")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[tissue %notin% c("limb","heart"), predActFromMTL_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_heart")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[tissue== "heart", predActFromMTL_VISTA_random := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_neuralTube")
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]
meta[tissue== "limb", predActFromMTL_VISTA_random := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_neuralTube")
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]
meta[tissue %notin% c("limb","heart"), predActFromMTL_VISTA_random := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_heart")
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]

# Melt into long format (each condition is linked to one file) ----
clean <- melt(meta,
              id.vars = c("ID", "tissue", "fold", "replicate"),
              value.name = "path")
clean[, c("file.type", "dataset", "set"):= tstrsplit(variable, "_")]
clean[, dataset:= ifelse(dataset=="VISTA", "activity", "accessibility")]
clean[, file.type:= paste0(file.type, "_file")]
clean[, file.type:= factor(file.type, unique(file.type))]
clean <- clean[file.exists(path)]
clean[, path:= normalizePath(path)]

# Melt into long format (each condition is linked to one file) ----
clean[, size:= file.size(path)/1e6]




# Dcast and save ----
final <- dcast(clean, ID+dataset+set+tissue+fold+replicate~file.type, value.var = "path")

# Save ----
saveRDS(final,
        "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/paper_metadata_v1.rds")
