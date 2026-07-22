setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
mot <- readRDS("Rdata/annotated_PWMs.rds")

# For each tissue ----
for(tiss in c("heart", "limb", "midbrain")){
  # Import sequences ----
  seq <- readRDS(paste0("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_", tiss, ".rds"))
  
  # For each background ----
  for(bg in c("genome", "subject")) {
    output.file <- paste0("db/motifs/motif_counts_", tiss, "_sequences_0.0001_", bg, "_bg.rds")
    if(!file.exists(output.file)) {
      # Compute motif counts ----
      counts <- vl_motifCounts(sequences = seq$sequence,
                               pwm_log_odds = mot$pwm,
                               genome = "mm10",
                               bg = bg,
                               p.cutoff = 1e-4)
      saveRDS(counts,
              output.file)
    }
  }
  
  # Selected labels ----
  sel <- c("vista_ts", "rdm", "vista_inactive_all_tissues",
           "evo.ini", "evo.act", "evo.act.acc",
           "ledidi_10_7_ini", "ledidi_10_9_ini", "ledidi_10_12_ini", "ledidi_10_14_ini", "ledidi_12_14_ini",
           "ledidi_10_7", "ledidi_10_9", "ledidi_10_12", "ledidi_10_14", "ledidi_12_14")
  
  # Compute motif enrichment versus random genomic sequences ----
  enr.file <- paste0("db/motifs/motif_enrich_", tiss, "_vs_rdm_seq_fisher.rds")
  if(!file.exists(enr.file)) {
    # Import counts
    counts <- readRDS(paste0("db/motifs/motif_counts_", tiss, "_sequences_0.0001_genome_bg.rds"))
    split.counts <- split(counts[seq$label %in% sel,], seq[label %in% sel, label])
    # Enrichment
    enr <- vl_motifEnrich(counts = split.counts[names(split.counts) != "rdm"],
                          control.counts = split.counts[["rdm"]],
                          names = mot$cluster)
    # Add contributions
    contrib <- mot[, c("motif", paste0("acc.contrib_", tiss), paste0("enh.contrib_", tiss)), with= FALSE]
    setnames(contrib, c("motif", "acc.contrib", "enh.contrib"))
    enr <- merge(enr,
                  contrib,
                  by= "motif",
                  all.x= T)
    # Save
    saveRDS(enr, enr.file)
  }
}