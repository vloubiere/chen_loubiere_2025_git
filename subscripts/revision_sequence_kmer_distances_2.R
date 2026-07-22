setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")
require(data.table)
require(Biostrings)

# Import initialization and designed enhancer sequences
heart <- readRDS("Rdata/final_designed_enhancer_sequences_heart.rds")
heart <- heart[id %in% c(311, 726, 834, 890, 845) & label == "ledidi_12_14"]

limb <- readRDS("Rdata/final_designed_enhancer_sequences_limb.rds")
limb <- limb[id %in% c(1104, 1121, 51, 103, 112) & label == "ledidi_12_14"]

CNS <- readRDS("Rdata/final_designed_enhancer_sequences_midbrain.rds")
CNS <- CNS[id %in% c(1008, 734, 543, 169, 5) & label == "ledidi_12_14"]

dat <- rbindlist(list(heart = heart, limb = limb, CNS = CNS), idcol = "tissue")
dat[, label := paste0(tissue, ".validated.enhancer")]

# Import vista sequences
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
vista[, c("seqnames", "start", "end", "strand", "name") := importBed(ifelse(genome == "hg38", coor_hg38, coor_mm10))]
vista <- vista[, .(genome, seqnames, start, end, strand = "+", heart, limb, midbrain, Ntissues = heart + limb + midbrain)]
vista[, sequence := getBSsequence(bed = .SD, genome = genome), genome]
vista <- vista[Ntissues <= 1]

vista[, label := fcase(
  heart == 1 & limb == 0 & midbrain == 0, "vista.heart",
  heart == 0 & limb == 1 & midbrain == 0, "vista.limb",
  heart == 0 & limb == 0 & midbrain == 1, "vista.CNS",
  Ntissues == 0, "vista.inactive"
)]

dat <- rbind(dat, vista, fill = TRUE)

# Keep only valid DNA
dat <- dat[!is.na(sequence)]
dat[, sequence := toupper(as.character(sequence))]
dat <- dat[grepl("^[ACGT]+$", sequence)]

# Canonical k-mer frequencies
k <- 4
dna <- DNAStringSet(dat$sequence)

kmer_mat <- oligonucleotideFrequency(dna, width = k)
kmer_names <- colnames(kmer_mat)

rc_names <- as.character(reverseComplement(DNAStringSet(kmer_names)))
canonical_names <- ifelse(kmer_names < rc_names, kmer_names, rc_names)

# Collapse kmer and reverse-complement columns
canon_levels <- unique(canonical_names)
canon_mat <- sapply(canon_levels, function(x) {
  cols <- which(canonical_names == x)
  rowSums(kmer_mat[, cols, drop = FALSE])
})

canon_mat <- as.matrix(canon_mat)
colnames(canon_mat) <- canon_levels
rownames(canon_mat) <- paste0(dat$label, "_", seq_len(nrow(dat)))

# Normalize to frequencies
mat <- canon_mat / rowSums(canon_mat)

# Remove zero columns and very rare kmers
keep1 <- colSums(mat) > 0
mat <- mat[, keep1, drop = FALSE]

keep2 <- colMeans(mat > 0) >= 0.05
mat <- mat[, keep2, drop = FALSE]

# Hellinger transform: often works better for compositional count/frequency data
mat_hel <- sqrt(mat)

# Cosine distance
cosine_dist <- function(x) {
  x <- as.matrix(x)
  nr <- sqrt(rowSums(x^2))
  sim <- tcrossprod(x) / outer(nr, nr)
  sim[is.na(sim)] <- 0
  as.dist(1 - sim)
}

d <- cosine_dist(mat_hel)

# PCoA
pcoa <- cmdscale(d, k = 2, eig = TRUE)
coord <- as.data.table(pcoa$points)
setnames(coord, c("V1", "V2"))

labels <- factor(
  dat$label,
  levels = c(
    "CNS.validated.enhancer",
    "heart.validated.enhancer",
    "limb.validated.enhancer",
    "vista.CNS",
    "vista.heart",
    "vista.limb",
    "vista.inactive"
  )
)

cols <- c("limegreen", "cornflowerblue", "tomato", "limegreen", "cornflowerblue", "tomato", "grey")[labels]

final <- data.table(
  V1 = coord$V1,
  V2 = coord$V2,
  label = labels,
  col = cols
)

final[, validated := grepl("validated", as.character(label))]

vl_par()
final[, plot(
  V1, V2,
  pch = 21,
  bg = adjustcolor(col, .5),
  col = ifelse(validated, "black", NA),
  xlab = "PCoA1",
  ylab = "PCoA2"
)]

# Optional hierarchical clustering
hc <- hclust(d, method = "ward.D2")
plot(hc)
