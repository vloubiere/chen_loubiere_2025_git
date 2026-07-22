setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# Import vista sequences ----
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
vista[, c("seqnames", "start", "end", "strand", "name"):= importBed(ifelse(genome=="hg38", coor_hg38, coor_mm10))]
vista <- vista[, .(genome, seqnames, start, end, strand= "+", heart, limb, midbrain, Ntissues= heart+limb+midbrain)]
vista <- rbindlist(
  list(
    heart= vista[heart==1],
    limb= vista[limb==1],
    CNS= vista[midbrain==1],
    Inactive= vista[heart+limb+midbrain==0] 
  ),
  idcol = "tissue"
)

# Design a set of randomly sampled genomic sequences for VISTA ----
set.seed(1)
rdm.vista <- vlite::randomRegionsBSgenome("mm10", widths = vista[, end-start+1])
rdm.vista[, genome:= "mm10"]
rdm.vista[, strand:= "+"]
rdm.vista$width <- NULL
rdm.vista[, tissue:= "inactive"]

# Import all designed enhancer sequences ----
heart <- readRDS("Rdata/final_designed_enhancer_sequences_heart.rds")
heart <- heart[label=="ledidi_12_14"]
limb <- readRDS("Rdata/final_designed_enhancer_sequences_limb.rds")
limb <- limb[label=="ledidi_12_14"]
CNS <- readRDS("Rdata/final_designed_enhancer_sequences_midbrain.rds")
CNS <- CNS[label=="ledidi_12_14"]
design <- rbindlist(list(heart= heart, limb= limb, CNS= CNS), idcol = "tissue")

# Import ATAC sequences ----
atac <- readRDS("db/peaks/bulkENCODE_confident_ATAC_peaks.rds")
atac[, tissue:= class]
atac[, genome:= "mm10"]
atac[, class:= "atac"]
atac[, strand:= "+"]
atac <- atac[tissue %in% c("heart", "midbrain", "limb")]

# Add randomly sampled genomic sequences of 1kb ----
set.seed(1)
rdm <- vlite::randomRegionsBSgenome("mm10", widths = rep(1001, 1200))
rdm[, genome:= "mm10"]
rdm[, strand:= "+"]
rdm$width <- NULL
rdm[, tissue:= "inactive"]

# Merge and extract sequences ----
dat <- rbindlist(
  list(
    vista= vista,
    rdm.vista= rdm.vista,
    design= design,
    atac= atac,
    rdm= rdm
  ),
  use.names= TRUE,
  fill= TRUE,
  idcol = "class"
)
dat[is.na(sequence), sequence:= getBSsequence(bed = .SD, genome= genome), genome]

# Count motifs ----
mot.db <- readRDS("../../motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")
mot <- vl_motifCounts(
  sequences = dat$sequence,
  pwm_log_odds = mot.db$pwms_log_odds,
  bg= "genome",
  genome = "mm10",
  p.cutoff = 1e-4
)

# Compute enrichment versus randomly sampled genomic sequences ----
dat[!grepl("rdm", class), .N, .(class, tissue)]
enr <- dat[!grepl("rdm", class), {
  ctl <- if(class=="vista")
    mot[dat$class=="rdm.vista",] else
      mot[dat$class=="rdm",]
  vl_motifEnrich(mot[dat$class==class & dat$tissue==tissue,], ctl)
}, .(class, tissue)]

# Sanity checks ----
enr[name=="GATA__GATA3_MA0037.3"]
enr[name=="HD/20__OTX2_MOUSE.H11MO.0.A"]
enr[name=="Ebox/CAGCTG__TWIST1_MA1123.1"]

# Save
saveRDS(
  enr,
  "db/motifs/revision_motif_enrich_atac_vista_designed_vs_rdm_genomic.rds"
)
