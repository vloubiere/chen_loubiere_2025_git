setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# PREPARE DATA FOR TRAINING ---------------------------------------------------------------------
# Functions
file.edit("git_deepATAC/function/augmentation_function_tiling_sliding_window.R") # For ATAC-Seq peaks
file.edit("git_deepATAC/function/compute_AUC.R") # Compute AUC
# Clean VISTA tiles
file.edit("git_deepATAC/subscripts/clean_VISTA_tiles.R")
# Define folds for ATAC/transfer learning with bulk ATAC-Seq dataset
file.edit("git_deepATAC/subscripts/confident_bulkATAC_peaks.R") # ATAC clustering
file.edit("git_deepATAC/subscripts/define_bulkATAC_control_regions.R") # Non-overlapping controls
file.edit("git_deepATAC/subscripts/nonOverlapping_bulkATAC_folds.R") # Split into non-overlapping train/valid/test
file.edit("git_deepATAC/subscripts/augment_bulkATAC_regions.R") # Augmentation
file.edit("git_deepATAC/subscripts/save_files_for_bulkATAC_training.R") # Save .fa .txt and .bed files
file.edit("git_deepATAC/subscripts/sanity_check.R") # Check that vista test set absent from train/valid sets
file.edit("git_deepATAC/subscripts/VISTA_chr18_PPV_bins.R") # Select chr18 bins for PPV sampling
# Make test set similar between approaches
file.edit("git_deepATAC/subscripts/select_comparable_test_sets.R") # Test set using only best design (to be able to compare between approaches)
# CNS model sets
file.edit("git_deepATAC/subscripts/CNS_combined_bw_track.R")

# PAPER ANALYSES --------------------------------------------------------------------------------
# Screenshot Figure 1 ----
file.edit("git_deepATAC/subscripts/screenshot_MEF2B_Fig1.R")

# Make clean metadata file ----
file.edit("git_deepATAC/subscripts/clean_metadata_paper.R") # Newly trained models (repdoducibility test, not used in paper)
file.edit("git_deepATAC/subscripts/clean_metadata_paper_v2.R") # Models used in paper (used for LEDIDI design)
file.edit("git_deepATAC/subscripts/sanity_check_overlap_training_test_validation_leakage.R") # Checks

# ACCESSIBILITY MODELS  -----------------------
# Compare 20 different approaches ATAC-Seq
file.edit("git_deepATAC/subscripts/heatmap_PCC_comparison_20_model_approaches.R")
file.edit("git_deepATAC/subscripts/heatmap_PPV_comparison_20_model_approaches.R")
# PCC scatterplot
file.edit("git_deepATAC/subscripts/scatterplot_PCC_accessibility_models_per_tissue.R")
file.edit("git_deepATAC/subscripts/scatterplot_PCC_accessibility_models_delta.R")
# Heatmap observed versus expected
file.edit("git_deepATAC/subscripts/heatmap_observed_vs_predicted_accessibility_models.R")
# Screenshot examples and PCC on chromosome 18
file.edit("git_deepATAC/subscripts/generate_bw_tracks_accessibility_predictions_chr18.R") # Compute tracks
file.edit("git_deepATAC/subscripts/observed_vs_predicted_PCC_chr18_peaks_union.R")
file.edit("git_deepATAC/subscripts/screenshot_example_accessbility_models.R")

# TRANSFER LEARNING ---------------------------
# PPV
file.edit("git_deepATAC/subscripts/PPV.R")
file.edit("git_deepATAC/subscripts/PPV_control_shuffled_models.R")
# Call motif positions
file.edit("git_deepATAC/subscripts/create_non_redundant_motif_set_Jeff.R")
file.edit("git_deepATAC/subscripts/collapsed_test_set_coor.R")
file.edit("git_deepATAC/subscripts/compute_motif_positions_collapsed_test_set_coor.R")
# Contribution scores motifs
file.edit("git_deepATAC/subscripts/VISTA_mean_motif_contrib_per_tissue.R")
file.edit("git_deepATAC/subscripts/ATAC_mean_motif_contrib_per_tissue.R")
file.edit("git_deepATAC/subscripts/compute_mean_contrib_per_motif.R")
# Motif contrib Heatmaps
file.edit("git_deepATAC/subscripts/heatmap_contrib_zscore_3_tissues.R")
file.edit("git_deepATAC/subscripts/heatmap_contrib_zscore_3_tissues_single_TF.R")
file.edit("git_deepATAC/subscripts/heatmap_contrib_zscore_all_tissues.R")
file.edit("git_deepATAC/subscripts/heatmap_contrib_wilcox.R")
# TFs GO enrichment
file.edit("git_deepATAC/subscripts/TFs_GO_enrichment.R")
# TFs single cell analysis
file.edit("git_deepATAC/subscripts/subset_single_cell_data.R")
file.edit("git_deepATAC/subscripts/AUCell_TF_clusters.R")
file.edit("git_deepATAC/subscripts/over_representation_TF_clusters.R")
