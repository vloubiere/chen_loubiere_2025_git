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

# Heatmap selected enhancers ----
file.edit("git_deepATAC/subscripts/compare_ledidi_parameters.R")
file.edit("git_deepATAC/subscripts/motifs_contribution_heart.R")
file.edit("git_deepATAC/subscripts/best_ledidi_parameters.R")
file.edit("git_deepATAC/subscripts/best_ledidi_parameters_heatmap.R")
file.edit("git_deepATAC/subscripts/known_heart_motifs.R")
file.edit("git_deepATAC/subscripts/all_known_motifs.R")

# Final selection enhancers ! ----
file.edit("git_deepATAC/subscripts/set_cutoffs_predicted_activity.R")
file.edit("git_deepATAC/subscripts/compute_motif_counts_and_enrichments.R")
file.edit("git_deepATAC/subscripts/compare_design_approaches.R")
file.edit("git_deepATAC/subscripts/heatmap_motif_counts_per_kb.R")
file.edit("git_deepATAC/subscripts/heatmap_compare_tissues.R")
file.edit("git_deepATAC/subscripts/best_ledidi_parameters_contrib.R")

# Final figures Evgeny ----
file.edit("git_deepATAC/subscripts/create_clean_list_sequences.R")
file.edit("git_deepATAC/subscripts/set_cutoffs_predicted_activity_evgeny.R")
file.edit("git_deepATAC/subscripts/annotate_PWMs.R")
file.edit("git_deepATAC/subscripts/compute_motif_counts_and_enrichments.R")
file.edit("git_deepATAC/subscripts/compare_design_approaches_evgeny.R")
file.edit("git_deepATAC/subscripts/heatmap_motif_counts_per_kb_evgeny.R")
file.edit("git_deepATAC/subscripts/heatmap_compare_tissues_evgeny.R")
file.edit("git_deepATAC/subscripts/cluster_filtered_sequences_evgeny.R")
file.edit("git_deepATAC/subscripts/contribution_tracks_selected_sequences_evgeny.R")

# PAPER FIRST VERSION ----
file.edit("git_deepATAC/subscripts/clean_metadata_paper.R")
file.edit("git_deepATAC/subscripts/PCC_paper.R")

