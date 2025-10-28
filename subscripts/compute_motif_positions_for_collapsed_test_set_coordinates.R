setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import collapsed test sets ----
bed <- importBed("db/bed/collapsed_test_set_paper.bed")

# Extract sequence ----
bed[, genome:= tstrsplit(name, "__", keep= 2)]
bed[, sequence:= getBSsequence(.SD, genome), genome]
seqs <- bed$sequence
names(seqs) <- bed$name

# Import motifs db files ----
mot <- readRDS("/groups/stark/vloubiere/motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")

# Compute motifs positions ----
mot.pos.file <- "db/motifs/non_redundant_motif_pos_test_set_paper_subject_0.0001.rds"
if(!file.exists(mot.pos.file)) {
  pos <- vl_motifPos(
    sequences = seqs,
    pwm_log_odds = mot$pwms_log_odds,
    bg = "subject",
    p.cutoff = 5e-05
  )
  # Remove empty sequences
  pos <- pos[!is.na(mot.count)]
  saveRDS(pos, mot.pos.file)
  # Format
  bed <- motifPosToBed(pos)
  bed[, c("seqnames", "genome"):= tstrsplit(seqlvls, "__")]
  bed[, c("seqnames", "region.start", "region.end"):= importBed(seqnames)[, .(seqnames, start, end)]]
  bed[, start:= region.start+start-1]
  bed[, end:= region.start+end-1]
  bed <- bed[, .(genome, seqnames, start, end, name= motif, score, strand)]
  # Collapse each motif coor
  coll <- bed[, collapseBed(.SD), .(genome, name)]
  saveRDS(coll,
          "db/motifs/non_redundant_motif_pos_test_set_paper_subject_0.0001_collapsed.rds")
}


