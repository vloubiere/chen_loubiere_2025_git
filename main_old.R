setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# OLD ----------------------------------------------------------------------------------------------------------------#
file.edit("git_deepATAC/subscripts/final_heatmap_selected_enhancers_old.R")
# ATAC test sets
file.edit("git_deepATAC/subscripts/TS_quantif_test_sets_ATAC_models.R") # tissueSpecific vs specificClosed
file.edit("git_deepATAC/subscripts/smoothScatter_test_sets_ATAC_models.R") # Correlation obs vs. mean pred
file.edit("git_deepATAC/subscripts/compare_ATAC_moldels_perf_test_set.R") # Correlation obs vs. mean pred
# ATAC chr18
file.edit("git_deepATAC/subscripts/assess_reproducibility_ATAC_chr18_and_output_mean_predictions.R") # PCC 
file.edit("git_deepATAC/subscripts/smoothScatter_chr18_ATAC_models.R") # Correlation obs vs. mean pred
# VISTA test sets and chr18
file.edit("git_deepATAC/subscripts/assess_reproducibility_VISTA_and_output_mean_predictions.R") # ROC AUC
file.edit("git_deepATAC/subscripts/mean_predictions_chr18_PPV_sampling.R") # ROC AUC
file.edit("git_deepATAC/subscripts/PPV_VISTA_models.R") # PPV mean 

# Old model 20240807 ----
file.edit("git_deepATAC/subscripts/retrieve_shenzhi_ATAC_models.R")
file.edit("git_deepATAC/subscripts/diagnostic_ATAC_models.R")
file.edit("git_deepATAC/subscripts/compare_PCC_between_folds_and_reps.R")
file.edit("git_deepATAC/subscripts/Heatmap_compare_observed_expected_zscores.R")
file.edit("git_deepATAC/subscripts/Scatterplot_observed_predicted_testSets.R")
file.edit("git_deepATAC/subscripts/Scatterplot_observed_predicted_chr11.R")
file.edit("git_deepATAC/subscripts/diagnostic_plots_models.R") # Not used

# Training using sn-ATAC-Seq ----
file.edit("git_deepATAC/subscripts/confident_sn_ATAC_Seq_peaks.R")

# Diagnostic VISTA models ----
file.edit("git_deepATAC/subscripts/retrieve_shenzhi_VISTA_models.R")
file.edit("git_deepATAC/subscripts/ROC_AUC_augmentation_VISTA.R")
file.edit("git_deepATAC/subscripts/CG_content_impact.R")
file.edit("git_deepATAC/subscripts/compare_old_new_VISTA_AUC.R") # Not used
file.edit("git_deepATAC/subscripts/check_augmentation_impact.R") # Not used

# Count enhancers vista database ----
file.edit("git_deepATAC/subscripts/barplot_VISTA_tiles_per_tissue.R")

# Unused ----
file.edit("git_deepATAC/subscripts/save_ATAC_sequences_and_coverage_for_each_tissue.R") # Prepare ATAC-Seq training
file.edit("git_deepATAC/subscripts/save_VISTA_sequences_and_activity_for_each_tissue.R") # Prepare VISTA transfer learning

