setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")
require(stringdist)

# Distance file
dist.file <- "db/sequence_distances/levenshtein_distances.rds"

if(!file.exists(dist.file)) {
  
  # Import initialization and designed enhancer sequences ----
  heart <- readRDS("Rdata/final_designed_enhancer_sequences_heart.rds")
  heart <- heart[(id %in% c(311, 726, 834, 890, 845) & label=="ledidi_12_14") | label=="ledidi_12_14_ini"]
  limb <- readRDS("Rdata/final_designed_enhancer_sequences_limb.rds")
  limb <- limb[(id %in% c(1104,  1121,  51,  103,  112) & label=="ledidi_12_14") | label=="ledidi_12_14_ini"]
  CNS <- readRDS("Rdata/final_designed_enhancer_sequences_midbrain.rds")
  CNS <- CNS[(id %in% c(1008, 734, 543, 169, 5) & label=="ledidi_12_14") | label=="ledidi_12_14_ini"]
  design <- rbindlist(list(heart= heart, limb= limb, CNS= CNS), idcol = "tissue")
  design[, label:= switch(
    label,
    "ledidi_12_14"= "validated.enhancer",
    "ledidi_12_14_ini"= "initialization.seq"
  ), label]
  design <- design[, .(tissue, label, init.id= id, sequence)]
  
  # Import vista sequences ----
  vista <- readRDS("db/peaks/vista_tiles_clean.rds")
  vista[, c("seqnames", "start", "end", "strand", "name"):= importBed(ifelse(genome=="hg38", coor_hg38, coor_mm10))]
  vista <- vista[, .(genome, seqnames, start, end, strand= "+", heart, limb, midbrain, Ntissues= heart+limb+midbrain)]
  vista <- rbind(vista, data.table::copy(vista)[, strand:= "-"])
  vista[, sequence:= getBSsequence(bed = .SD, genome= genome), genome]
  vista <- list(
    heart= vista[, .(label= ifelse(heart==1, "vista.enhancer", "vista.inactive"), sequence)],
    limb= vista[, .(label= ifelse(limb==1, "vista.enhancer", "vista.inactive"), sequence)],
    CNS= vista[, .(label= ifelse(midbrain==1, "vista.enhancer", "vista.inactive"), sequence)]
  )
  vista <- rbindlist(vista, idcol = "tissue")
  
  # Add random genomic sequences
  set.seed(1)
  rdm <- randomRegionsBSgenome("mm10", rep(1001, 1000))
  rdm[, label:= "rdm.genomic.seq"]
  rdm[, sequence:= vlite::getBSsequence(rdm, genome = "mm10")]
  rdm <- lapply(1:3, function(x) rdm)
  names(rdm) <- c("heart", "limb", "CNS")
  rdm <- rbindlist(rdm, idcol = "tissue")
  
  # All pairwise combinations ----
  dat <- rbind(design, vista, rdm, fill= T)
  dat[, id:= .I]
  comb <- CJ(dat[label=="validated.enhancer", id], dat[, id])
  a <- dat[comb$V1]
  b <- dat[comb$V2]
  setnames(a, function(x) paste0(x, ".a"))
  setnames(b, function(x) paste0(x, ".b"))
  
  # Select relevant combinations ----
  comb <- cbind(a, b)
  comb <- comb[(id.a != id.b) & (tissue.a == tissue.b)]
  comb[label.b=="initialization.seq" & init.id.a!=init.id.b, label.b:= "rdm.init.seq"]
  unique(comb[, .(label.a, label.b)])
  comb[, label.b:= factor(label.b, c("initialization.seq", "rdm.init.seq", "rdm.genomic.seq", "vista.inactive", "vista.enhancer", "validated.enhancer"))]
  setorderv(comb, "label.b")
  unique(comb[, .(label.a, label.b)])
  
  # Compute distances ----
  comb[, dist:= stringdist(sequence.a, sequence.b, method= "lv")]
  comb[, nchar.a:= nchar(sequence.a)]
  comb[, nchar.b:= nchar(sequence.b)]
  comb[, nchar.max:= apply(.SD, 1, max), .SDcols= c("nchar.a", "nchar.b")]
  comb[, norm:= dist/nchar.max]
  
  # Save ----
  saveRDS(comb, dist.file)
} else
  comb <- readRDS(dist.file)

# Only retain smaller distance for each sequence ----
dist <- comb[, .SD[which.min(dist)], .(tissue.a, sequence.a, label.b)]

# Plot ----
pdf("pdf/_revision/barplot_lv_distance_validated_enh.pdf", 4.5, 2.25)
vl_par(mai= c(.9, .9, .5, .5))
ls <- split(dist$norm, dist[, .(label.b, tissue.a)])
vl_barplot(
  ls,
  col= rep(adjustcolor(c("cornflowerblue", "limegreen", "tomato"), .3), each= nlevels(dist$label.b)),
  ylab= "Min. norm. lv distance"
)
dev.off()
