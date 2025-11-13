setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Functions -------------------------------------------------------------------------------------
file.edit("chen_loubiere_2025_git/function/augmentation_function_tiling_sliding_window.R") # For ATAC-Seq peaks

# PREPARE DATA FOR TRAINING ---------------------------------------------------------------------
# Clean VISTA tiles
file.edit("chen_loubiere_2025_git/subscripts/clean_VISTA_tiles.R")
# Cluster and annotate ATAC-Seq peaks !!! peak files have changed !!!
file.edit("chen_loubiere_2025_git/subscripts/clustering_confident_bulkATAC_peaks.R")
# Define inaccessible control regions (no overlap with ATAC-seq or VISTA) !!! peak files have changed !!!
file.edit("chen_loubiere_2025_git/subscripts/define_bulkATAC_control_regions.R")

# Define non-overlapping train/valid/test sets (one per fold) -----------------------------------
file.edit("chen_loubiere_2025_git/subscripts/nonOverlapping_bulkATAC_folds.R")

# DATA AUGMENTATION -----------------------------------------------------------------------------
# Augmentation
file.edit("chen_loubiere_2025_git/subscripts/augment_bulkATAC_regions.R") 
# Save .fa .txt and .bed files
file.edit("chen_loubiere_2025_git/subscripts/save_files_for_bulkATAC_training.R") 

# EXTRA TEST SETS FOR CV ------------------------------------------------------------------------
# Select chr18 bins for PPV sampling
file.edit("chen_loubiere_2025_git/subscripts/VISTA_chr18_PPV_bins.R")
# Test set using only best design (to compare between approaches)
file.edit("chen_loubiere_2025_git/subscripts/select_comparable_test_sets.R")

# SELECTED SEQUENCES FOR VALIDATIONS (round #1) ------------------------------------------------
# Clean list of designed sequences (different EVO/LEDIDI designs...)
file.edit("chen_loubiere_2025_git/subscripts/create_clean_list_sequences.R") 
# BLAST/specificity filters for LEDIDI_12_14...
file.edit("chen_loubiere_2025_git/subscripts/set_cutoffs_predicted_activity_evgeny.R")

# Make clean metadata file ---------------------------------------------------------------------
file.edit("chen_loubiere_2025_git/subscripts/clean_metadata_paper_v3.R") # Models used in paper (used for LEDIDI design)

# PAPER FIGURES --------------------------------------------------------------------------------
# FIGURE 1 ----
# MEF2B Screenshot (Fig. 1a)
file.edit("chen_loubiere_2025_git/subscripts/screenshot_MEF2B_Fig1.R")
# Screenshot examples and PCC on chromosome 18 (Fig. 1b)
file.edit("chen_loubiere_2025_git/subscripts/generate_bw_tracks_accessibility_predictions_chr18.R") # Compute tracks
file.edit("chen_loubiere_2025_git/subscripts/observed_vs_predicted_PCC_chr18_peaks_union.R")
file.edit("chen_loubiere_2025_git/subscripts/screenshot_example_accessbility_models.R")
# Heatmap observed versus expected (Fig. 1c and S1a)
file.edit("chen_loubiere_2025_git/subscripts/heatmap_observed_vs_predicted_accessibility_models.R")
# PPV (Fig. 1d)
file.edit("chen_loubiere_2025_git/subscripts/PPV.R")

# FIGURE S1 ----
# Barplot % TSS overlaps ATAC-Seq peaks cluster (Fig. S1a, see Fig. 1c)
# PCC scatterplots accessibility models (S1b-c)
file.edit("chen_loubiere_2025_git/subscripts/scatterplot_PCC_accessibility_models_per_tissue.R")
file.edit("chen_loubiere_2025_git/subscripts/scatterplot_PCC_accessibility_models_delta.R")
# Overlap VISTA tiles (S1d)
file.edit("chen_loubiere_2025_git/subscripts/heatmap_overlap_VISTA_tiles.R")
# Overlap midbrain enhancers CNS sub regions (S1e)
file.edit("chen_loubiere_2025_git/subscripts/pie_chart_overlap_midbrain_sub_CNS_regions.R")
# PPV other CNS tissues midbrain model (S1f)
file.edit("chen_loubiere_2025_git/subscripts/PPV_midbrain_other_CNS_tissues.R")
# PPV other CNS tissues midbrain model (S1g)
file.edit("chen_loubiere_2025_git/subscripts/boxplot_predict_activity_midbrain_model_VISTA_enhancers.R")
# Negative sequences rejection (S1h)
file.edit("chen_loubiere_2025_git/subscripts/table_FALSE_POSITIVE_sequences.R")
# PPV acc. model and random initiation models (S1i)
file.edit("chen_loubiere_2025_git/subscripts/PPV_control_shuffled_models.R")
# PPV difference between TL models and others (S1j)
file.edit("chen_loubiere_2025_git/subscripts/PPV_control_shuffled_models_diff.R")

# FIGURE 2 ----
# Motif contrib Heatmap (Fig. 2a) ----
# Call motif positions
file.edit("chen_loubiere_2025_git/subscripts/create_non_redundant_motif_set_Jeff.R")
file.edit("chen_loubiere_2025_git/subscripts/collapsed_test_set_coordinates.R")
file.edit("chen_loubiere_2025_git/subscripts/compute_motif_positions_for_collapsed_test_set_coordinates.R")
# Contribution mean contribution score per motif instance
file.edit("chen_loubiere_2025_git/subscripts/VISTA_mean_motif_contrib_per_tissue.R")
file.edit("chen_loubiere_2025_git/subscripts/ATAC_mean_motif_contrib_per_tissue.R")
file.edit("chen_loubiere_2025_git/subscripts/compute_mean_contrib_per_motif.R")
# Motif contrib Heatmaps
file.edit("chen_loubiere_2025_git/subscripts/heatmap_contrib_zscore_3_tissues.R")
# Motif contribution tracks (Fig. 2b-d) ----
file.edit("chen_loubiere_2025_git/subscripts/contribution_tracks_designed_enh.R")
file.edit("chen_loubiere_2025_git/subscripts/contribution_tracks_designed_enh_new.R")

# FIGURE S2 ----
# TFs GO enrichment (S2a)
file.edit("chen_loubiere_2025_git/subscripts/TFs_GO_enrichment.R")
# TFs single cell analysis (S2b)
file.edit("chen_loubiere_2025_git/subscripts/subset_single_cell_data.R")
file.edit("chen_loubiere_2025_git/subscripts/AUCell_TF_clusters.R")
file.edit("chen_loubiere_2025_git/subscripts/over_representation_TF_clusters.R")
# Compare predicted activity between tissues
file.edit("chen_loubiere_2025_git/subscripts/compare_predicted_act_designed_seq_boxplot.R")
file.edit("chen_loubiere_2025_git/subscripts/barplot_percentage_specific_designed_seq.R")