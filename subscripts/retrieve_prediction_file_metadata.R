setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import metadata ----
meta <- readRDS("Rdata/metadata_ATACSeq_models_processed.rds")

# Retrieve ATAC predictions ----
ATACfolder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/model/ATAC_model/"
meta[, pred_ATAC_test:= file.path(ATACfolder, ID, tissue, paste0("results_", fold, "_", tissue, "_DeepSTARR2_", replicate), paste0(fold, "_sequences_test.fa_predictions_enhancer_Model.txt")), .(ID, tissue, fold, replicate)]
meta[, pred_ATAC_chr18:= file.path(ATACfolder, ID, tissue, paste0("results_", fold, "_", tissue, "_DeepSTARR2_", replicate), "whole_chr18.fa_predictions_enhancer_Model.txt"), .(ID, tissue, fold, replicate)]

# Retrieve VISTA predictions ----
VISTAfolder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/model/VISTA_model/"
meta[, pred_VISTA_test:= file.path(VISTAfolder, ID, tissue, paste0("results_", fold, "_", tissue, "_DeepSTARR2_", replicate), paste0(fold, "_sequences_test.fa_predictions_enhancer_Model.txt")), .(ID, tissue, fold, replicate)]
meta[, pred_VISTA_chr18:= file.path(VISTAfolder, ID, tissue, paste0("results_", fold, "_", tissue, "_DeepSTARR2_", replicate), "PPV_bins_chr18.fa_predictions_enhancer_Model.txt"), .(ID, tissue, fold, replicate)]
meta[, pred_VISTA_rdm:= file.path(VISTAfolder, ID, tissue, paste0("results_", fold, "_", tissue, "_DeepSTARR2_", replicate), "random_sequences_600k.fasta_predictions_enhancer_Model.txt"), .(ID, tissue, fold, replicate)]

# Save
saveRDS(meta,
        "Rdata/metadata_ATACSeq_model_predictions.rds")