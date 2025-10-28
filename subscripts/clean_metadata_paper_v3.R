# Set workpath ----
setwd("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
options(datatable.prettyprint.char = 20)

# Import models metadata ----
meta <- readxl::read_xlsx("clean_version/Rdata/metadata_ATACSeq_models.xlsx", skip = 4)
meta <- as.data.table(meta)

# For each approach, add tissues, folds and replicates ----
meta <- meta[, .(tissue= c("heart", "forebrain", "hindbrain", "midbrain", "neuralTube", "limb","CNS")), (meta)]
meta <- meta[, .(fold= paste0("fold0", 1:3)), (meta)]
meta <- meta[, .(replicate= paste0("rep", 1:2)), (meta)]
meta <- meta[, .(set= c("training", "validation", "test")), (meta)]

# Retrieve input files ----
# bed
meta[tissue=="CNS", bed_VISTA:= paste0("db/training_dataset/VISTA_20240901/", tissue, "/", fold, "_sequences_", set,".bed")]
meta[tissue!="CNS", bed_VISTA:= paste0("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/bed/",
                                       dataset,"/VISTA/", tissue, "/", fold, "_sequences_", set,".bed")]
meta[tissue=="CNS", bed_ATAC:= paste0("db/training_dataset/", Batch, "_bulkATAC_", augmentation, "Aug_", balancing, "Bal/",  tissue, "/", fold, "_sequences_", set,".bed")]
meta[tissue!="CNS", bed_ATAC:= paste0("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/bed/",
                                      dataset,"/ATAC/",  tissue, "/", augmentation, "/", balancing, "/", fold, "_sequences_", set,".bed")]
# fasta
meta[, fa_VISTA:= paste0("db/training_dataset/VISTA_20240901/",
                         tissue, "/", fold, "_sequences_", set,".fa")]
meta[, fa_ATAC:= paste0("db/training_dataset/",
                        Batch, "_bulkATAC_", augmentation, "Aug_", balancing, "Bal/", tissue, "/", fold, "_sequences_", set,".fa")]
# Observed (coverage for ATAC, 0/1 labels for VISTA)
meta[, obs_VISTA:= paste0("db/training_dataset/VISTA_20240901/",
                          tissue, "/", fold, "_sequences_activity_", set,".txt")]
meta[, obs_ATAC:= paste0("db/training_dataset/",
                         Batch, "_bulkATAC_", augmentation, "Aug_", balancing, "Bal/", tissue, "/", fold, "_sequences_activity_", set,".txt")]

# Dcast to have one set per column ----
meta <- dcast(meta,
              augmentation + balancing + weight + ID + tissue + fold + replicate ~ set,
              value.var = list("bed_VISTA", "bed_ATAC", "fa_VISTA", "fa_ATAC", "obs_VISTA", "obs_ATAC"))

# Retrieve files for the test set generated with the best design, to allow comparison between 20 approaches ---
meta[, obs_access_testBestDesign := paste0("db/training_dataset/model1_bulkATAC_",
                                           augmentation, "Aug_", balancing, "Bal/", tissue, "/", fold, "_sequences_testBestDesign.txt")]
meta[, fa_access_testBestDesign := paste0("db/training_dataset/model1_bulkATAC_",
                                          augmentation, "Aug_", balancing, "Bal/", tissue, "/", fold, "_sequences_testBestDesign.fa")]
meta[tissue!="CNS",
     bed_access_testBestDesign := paste0("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/bed/bulkATAC/ATAC/", 
                                         tissue, "/", augmentation, "/", balancing, "/", fold, "_sequences_testBestDesign.bed")]
meta[tissue=="CNS",
     bed_access_testBestDesign := paste0("db/training_dataset/model1_bulkATAC_",
                                         augmentation, "Aug_", balancing, "Bal/", tissue, "/", fold, "_sequences_testBestDesign.bed")]

# Retrieve files for the full chr18, used for screenshots and as an alternative test set ---
meta[, obs_access_chr18 := paste0("db/testing_dataset/chr18/", tissue, "_chr18_bins_access_obs.txt")]
meta[, fa_access_chr18 := paste0("db/testing_dataset/chr18/chr18_bins.fa")]
meta[, bed_access_chr18 := paste0("db/testing_dataset/chr18/", tissue, "_chr18_bins_access_obs.bed")]

# Add fasta files containing random sequences (used as extra test set for PPV) ----
meta[, fa_access_random:= "db/fasta/testing_dataset/random_sequences_600k.fasta"] # These are actually 300k!
meta[, fa_VISTA_random:= "db/fasta/testing_dataset/random_sequences_600k.fasta"] # These are actually 300k!

# Add fasta files containing negative control sequences (used as extra test set) ----
meta[, fa_access_NegGenomicRegions:= "db/fasta/testing_dataset/PPV_bins_chr18.fa"] # These are actually 300k!
meta[, fa_VISTA_NegGenomicRegions:= "db/fasta/testing_dataset/PPV_bins_chr18.fa"] # These are actually 300k!

# Add missing VISTA predictions ----
meta[, c("pred_VISTA_test", "pred_VISTA_random", "pred_VISTA_NegGenomicRegions"):= {
  # folder 
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  .(
    file.path(folder, paste0(fold, "_sequences_test.fa_predictions_enhancer_Model.txt")),
    file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt"),
    file.path(folder, "PPV_bins_chr18.fa_predictions_enhancer_Model.txt")
  )
}, .(ID, tissue, fold, replicate)]

# Retrieve accessibility prediction files (ATAC) ----
meta[, c("pred_access_test", "pred_access_random", "pred_access_chr18", "pred_access_testBestDesign"):= {
  # folder 
  folder <- paste0("result/model/ATAC_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  .(
    file.path(folder, paste0(fold, "_sequences_test.fa_predictions_enhancer_Model.txt")),
    file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt"),
    file.path(folder, "whole_chr18.fa_predictions_enhancer_Model.txt"),
    file.path(folder, paste0(fold, "_sequences_testBestDesign.fa_predictions_enhancer_Model.txt"))
  )
}, .(ID, tissue, fold, replicate)]

# Add contribution files ----
# Containing full test set (very large)
meta[, contrib_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  file.path(folder, paste0("/Model_", fold,"_sequences_test.fa_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5"))
}]
meta[, contrib_access_test := {
  folder <- paste0("result/model/ATAC_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  file.path(folder, paste0("/Model_", fold ,"_sequences_test.fa_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5"))
}]
# For ATAC, only use center bins (Note that there is no fasta file here)
meta[, bed_access_testCenterActBins:= paste0("db/training_dataset/model1_bulkATAC_",
                                             augmentation, "Aug_", balancing, "Bal/", tissue, "/", fold, "_sequences_testContrib.bed"), .(tissue, fold)]
meta[, fa_access_testCenterActBins:= paste0("db/training_dataset/model1_bulkATAC_",
                                            augmentation, "Aug_", balancing, "Bal/", tissue, "/", fold, "_sequences_testContrib.fa"), .(tissue, fold)]
meta[, contrib_access_testCenterActBins := {
  folder <- paste0("result/model/ATAC_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  file.path(folder, paste0("/Model_", fold ,"_sequences_testContrib.fa_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5"))
}]

# Add activity prediction from accessibility model (control) ----
meta[, predActFromAcc_VISTA_test := {
  folder <- paste0("result/model/ATAC_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "/act")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[, predActFromAcc_VISTA_random := {
  folder <- paste0("result/model/ATAC_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate)
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]

# Add activity prediction from randomly initialized models (control)  ----
meta[, predActFromRandomIni_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_init_random")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[, predActFromRandomIni_VISTA_random := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_init_random")
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]

# Add activity prediction from mismatched tissue TL models (controls) ----
meta[tissue %in% c("limb","heart"), predActFromMTL_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_midbrain")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[tissue %notin% c("limb","heart"), predActFromMTL_VISTA_test := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_heart")
  file.path(folder, paste0(fold,"_sequences_test.fa_predictions_enhancer_Model.txt"))
}]
meta[tissue %in% c("limb","heart"), predActFromMTL_VISTA_random := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_midbrain")
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]
meta[tissue %notin% c("limb","heart"), predActFromMTL_VISTA_random := {
  folder <- paste0("result/model/VISTA_model/", ID, "/", tissue, "/results_", fold, "_", tissue, "_DeepSTARR2_", replicate, "_transfer_heart")
  file.path(folder, "random_sequences_600k.fasta_predictions_enhancer_Model.txt")
}]

# Clean and melt into long format (each condition is linked to one file) ----
meta$augmentation <- meta$balancing <- meta$weight <- NULL
clean <- melt(meta,
              id.vars = c("ID", "tissue", "fold", "replicate"),
              value.name = "path")

# Retrieve info columns ----
clean[, c("file.type", "dataset", "set"):= tstrsplit(variable, "_")]
clean[, dataset:= ifelse(dataset=="VISTA", "activity", "accessibility")]
clean[, file.type:= paste0(file.type, "_file")]
clean[, file.type:= factor(file.type, unique(file.type))]
clean <- clean[file.exists(path)]
clean[, path:= normalizePath(path)]

# Dcast and SAVE ----
final <- dcast(clean, ID+dataset+set+tissue+fold+replicate~file.type, value.var = "path")
saveRDS(final,
        "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/paper_metadata_v3.rds")
