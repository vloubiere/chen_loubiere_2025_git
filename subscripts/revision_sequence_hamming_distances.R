setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")
require(stringdist)

# Distance file
dist.file <- "db/sequence_distances/hamming_distances.rds"

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
  # Tile
  vista[end-start+1<1001, c("start", "end"):= .(round(rowMeans(.SD))-500, round(rowMeans(.SD))+500), .SDcols= c("start", "end")]
  vista <- binBed(vista, bins.width = 1001, steps.width = 1)
  # Extract sequences
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
  
  # Compute distances ----
  dat <- rbind(design, vista, rdm, fill= T)
  dist <- dat[label=="validated.enhancer", {
    .s <- sequence
    .t <- tissue
    .c <- dat[tissue==.t, {
      .(init.id.b= init.id, dist= stringdist(.s, sequence, method= "hamming"))
    }, label]
    .c
  }, .(tissue, init.id.a= init.id, sequence)]
  # Remove similar sequences
  dist <- dist[!(label=="validated.enhancer" & init.id.a==init.id.b)]
  # Order labels
  dist[label=="initialization.seq" & init.id.a!=init.id.b, label:= "rdm.init.seq"]
  dist[, label:= factor(label, c("initialization.seq", "rdm.init.seq", "rdm.genomic.seq", "vista.inactive", "vista.enhancer", "validated.enhancer"))]
  
  # Only retain smaller distance for each sequence ----
  res <- dist[, .SD[which.min(dist)], .(tissue, label, sequence)]
  setorderv(res, "label")
  
  # Save ----
  saveRDS(res, dist.file)
} else
  dist <- readRDS(dist.file)

# Plot ----
pdf("pdf/_revision/barplot_hamming_distance_validated_enh.pdf", 4.5, 2.25)
vl_par(mai= c(.9, .9, .5, .5))
ls <- split(dist$dist, dist[, .(label, tissue)])
vl_barplot(
  ls,
  col= rep(adjustcolor(c("cornflowerblue", "limegreen", "tomato"), .3), each= nlevels(dist$label)),
  ylab= "Min. hamming distance"
)
dev.off()
