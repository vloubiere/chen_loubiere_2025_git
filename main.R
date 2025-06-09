setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Functions ----
file.edit("git_deepATAC/function/augmentation_function_tiling_sliding_window.R") # For ATAC-Seq peaks
file.edit("git_deepATAC/function/compute_AUC.R") # Compute AUC

# Clean VISTA tiles ----
file.edit("git_deepATAC/subscripts/clean_VISTA_tiles.R")

# Define folds for ATAC/transfer learning with bulk ATAC-Seq dataset ----
file.edit("git_deepATAC/subscripts/confident_bulkATAC_peaks.R") # ATAC clustering
file.edit("git_deepATAC/subscripts/define_bulkATAC_control_regions.R") # Non-overlapping controls
file.edit("git_deepATAC/subscripts/nonOverlapping_bulkATAC_folds.R") # Split into non-overlapping train/valid/test
file.edit("git_deepATAC/subscripts/augment_bulkATAC_regions.R") # Augmentation
file.edit("git_deepATAC/subscripts/save_files_for_bulkATAC_training.R") # Save .fa .txt and .bed files
file.edit("git_deepATAC/subscripts/sanity_check.R") # Check that vista test set absent from train/valid sets
file.edit("git_deepATAC/subscripts/VISTA_chr18_PPV_bins.R") # Select chr18 bins for PPV sampling

# Diagnostics ATAC/VISTA models gradient ----
file.edit("git_deepATAC/subscripts/retrieve_prediction_file_metadata.R") # Retrieve files
file.edit("git_deepATAC/subscripts/model1_results_table.R") # Make full results table, not used
# Compute performance metrics
file.edit("git_deepATAC/subscripts/compute_ATAC_PCC_per_replicate_and_mean_predictions.R") # ATAC performance
file.edit("git_deepATAC/subscripts/compute_VISTA_AUC_mPCC_PPV_per_replicate_and_mean_predictions.R") # VISTA performance
# Plot metrics per replicate
file.edit("git_deepATAC/subscripts/plot_metrics_per_replicate.R")
# ATAC-Seq performance
file.edit("git_deepATAC/subscripts/smoothScatter_mean_ATAC_predictions.R") # ATAC-Seq obs/exp
# VISTA performance
file.edit("git_deepATAC/subscripts/density_VISTA_predictions.R") # Predictions density
file.edit("git_deepATAC/subscripts/PPV_VISTA_predictions.R") # Postive Predicted Values analysis
file.edit("git_deepATAC/subscripts/TPR_VISTA_predictions.R") # TRUE positive rate
pdftools::pdf_combine(c("pdf/model1_performance_metrics_per_replicate.pdf",
                        "pdf/model1_ATAC_smoothScatter_mean_predictions.pdf",
                        "pdf/model1_VISTA_predictions_density_plots.pdf",
                        "pdf/model1_VISTA_PPV_curves.pdf",
                        "pdf/model1_VISTA_TPR_tile_counts.pdf"),
                      output = "pdf/model1_performance_merged.pdf")

# ATAC-Seq prediction snATAC-Seq ----
file.edit("git_deepATAC/subscripts/ATAC_seq_predictive_of_enhancer_act.R")

# Motif enrichment vista enhancers -----
file.edit("git_deepATAC/subscripts/motif_enrichment_VISTA_enhancers.R")

# Count motifs in designed sequences -----
file.edit("git_deepATAC/subscripts/count_motifs_designed_sequences.R")

# Clustering of motifs in designed enhancers vs. vista enhancers -----
file.edit("git_deepATAC/subscripts/clustering_motifs_designed_sequences.R")

# Call seqlets from VISTA enhancers ----
file.edit("git_deepATAC/subscripts/call_seqlets_heart_VISTA_enhancers.R")



# OLD ----------------------------------------------------------------------------------------------------------------#
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

