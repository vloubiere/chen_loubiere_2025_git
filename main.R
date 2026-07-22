setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Functions -------------------------------------------------------------------------------------
file.edit("chen_loubiere_2025_git/function/augmentation_function_tiling_sliding_window.R") # For ATAC-Seq peaks
file.edit("chen_loubiere_2025_git/function/compute_AUC.R") # Compute AUC

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
# Check that vista test set absent from train/valid sets
file.edit("chen_loubiere_2025_git/subscripts/sanity_check.R") 

# EXTRA TEST SETS FOR CV ------------------------------------------------------------------------
# Select chr18 bins for PPV sampling
file.edit("chen_loubiere_2025_git/subscripts/VISTA_chr18_PPV_bins.R")
# Test set using only best design (to compare between approaches)
file.edit("chen_loubiere_2025_git/subscripts/select_comparable_test_sets.R")
# CNS model sets
file.edit("chen_loubiere_2025_git/subscripts/CNS_combined_bw_track.R")

# SELECTED SEQUENCES FOR VALIDATIONS (round #1) ------------------------------------------------
# Clean list of designed sequences (different EVO/LEDIDI designs...)
file.edit("chen_loubiere_2025_git/subscripts/create_clean_list_sequences.R") 
# BLAST/specificity filters for LEDIDI_12_14...
file.edit("chen_loubiere_2025_git/subscripts/set_cutoffs_predicted_activity_evgeny.R")

# DESIGN ENHANCERS SPECIFIC SUB-CNS regions ----------------------------------------------------
file.edit("chen_loubiere_2025_git/subscripts/compare_act_acc_LEDIDI_settings_across_tissues.R")
file.edit("chen_loubiere_2025_git/subscripts/SOX3_impact.R")

# Make clean metadata file ---------------------------------------------------------------------
file.edit("chen_loubiere_2025_git/subscripts/clean_metadata_paper_v3.R") # Models used in paper (used for LEDIDI design)
file.edit("chen_loubiere_2025_git/subscripts/sanity_check_overlap_training_test_validation_leakage.R") # Sanity check

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
file.edit("chen_loubiere_2025_git/subscripts/barplot_overlap_VISTA_tiles.R")
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

# DELETED FIGURES -----------------------
# Compare 20 different approaches
file.edit("chen_loubiere_2025_git/subscripts/heatmap_PCC_comparison_20_model_approaches.R")
file.edit("chen_loubiere_2025_git/subscripts/heatmap_PPV_comparison_20_model_approaches.R")

# REVISION Nature Genetics ----------------------------------------------------------------------

# Transfer learning using chromatin-inferred enhancer labels ----
# ChIP-seq peaks
file.edit("chen_loubiere_2025_git/subscripts/ChIPseq_peaks_rds.R")
file.edit("chen_loubiere_2025_git/subscripts/ChIPseq_QC_check_overlaps_between_peaks.R") # Sanity check
# Define labels for transfer-learning
file.edit("chen_loubiere_2025_git/subscripts/define_VISTA_labels_based_on_enhancers_hallmarks.R") # VISTA (TL1 -> TL7, not used!)
file.edit("chen_loubiere_2025_git/subscripts/define_ATAC_labels_based_on_enhancers_hallmarks.R") # ATAC peaks (TL8 -> TL14, not used!)
file.edit("chen_loubiere_2025_git/subscripts/enhancer_hallmarks_CV_on_VISTA.R")
file.edit("chen_loubiere_2025_git/subscripts/enhancer_hallmarks_number_of_sequences.R")

# Size distribution vista enhancers ----
file.edit("chen_loubiere_2025_git/subscripts/revision_vista_size_distrib.R")

# Sequence distance ----
# Levenshtein
file.edit("chen_loubiere_2025_git/subscripts/revision_sequence_levenshtein_distances.R")
# Hamming
file.edit("chen_loubiere_2025_git/subscripts/revision_sequence_hamming_distances.R")
# Motifs
file.edit("chen_loubiere_2025_git/subscripts/motif_enrichment_vista_and_designed_sequences.R")
file.edit("chen_loubiere_2025_git/subscripts/count_motifs_validated_enhancers.R")
file.edit("chen_loubiere_2025_git/subscripts/revision_design_all_motifs.R")
file.edit("chen_loubiere_2025_git/subscripts/revision_design_top_motifs.R")
# K-mer (not used)
file.edit("chen_loubiere_2025_git/subscripts/revision_sequence_kmer_distances.R") # K-mers not used for now
file.edit("chen_loubiere_2025_git/subscripts/revision_sequence_kmer_distances_2.R") # K-mers not used for now

# TL reweighting reproducibility ----
file.edit("chen_loubiere_2025_git/subscripts/boxplot_motif_contrib_reproducibility.R")

# PCC chr18 vs other test set ----
file.edit("chen_loubiere_2025_git/subscripts/scatterplot_PCC_accessibility_models_chr18.R")
file.edit("chen_loubiere_2025_git/subscripts/scatterplot_PCC_accessibility_models_other_chromosomes.R")

# Activity initilization sequences ----
file.edit("chen_loubiere_2025_git/subscripts/predicted_activity_init_sequences.R")
file.edit("chen_loubiere_2025_git/subscripts/predicted_activity_random_genomic_sequences.R")

# Exploratory analyses (not used) ---------------------------------------------------------------
# TRUE positive rate within VISTA enhancers using  ----
file.edit("chen_loubiere_2025_git/subscripts/tissue_specific_ATAC_clusters.R")
file.edit("chen_loubiere_2025_git/subscripts/tissue_specific_ATAC_peaks.R")
file.edit("chen_loubiere_2025_git/subscripts/tissue_specific_HTM_peaks.R")
file.edit("chen_loubiere_2025_git/subscripts/classificaiton_tissue_specific_RNAs.R")
file.edit("chen_loubiere_2025_git/subscripts/tissue_specific_RNA_genes.R")

