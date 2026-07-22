setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")
require(stringdist)

# Import initialization and designed enhancer sequences
heart <- readRDS("Rdata/final_designed_enhancer_sequences_heart.rds")
heart <- heart[id %in% c(311, 726, 834, 890, 845) & label=="ledidi_12_14"]
limb <- readRDS("Rdata/final_designed_enhancer_sequences_limb.rds")
limb <- limb[id %in% c(1104,  1121,  51,  103,  112) & label=="ledidi_12_14"]
CNS <- readRDS("Rdata/final_designed_enhancer_sequences_midbrain.rds")
CNS <- CNS[id %in% c(1008, 734, 543, 169, 5) & label=="ledidi_12_14"]
dat <- rbindlist(list(heart= heart, limb= limb, CNS= CNS), idcol = "tissue")
dat[, label:= paste0(tissue, ".validated.enhancer")]

# Import vista sequences
vista <- readRDS("db/peaks/vista_tiles_clean.rds")
vista[, c("seqnames", "start", "end", "strand", "name"):= importBed(ifelse(genome=="hg38", coor_hg38, coor_mm10))]
vista <- vista[, .(genome, seqnames, start, end, strand= "+", heart, limb, midbrain, Ntissues= heart+limb+midbrain)]
vista <- rbind(vista, data.table::copy(vista)[, strand:= "-"])
vista[, sequence:= getBSsequence(bed = .SD, genome= genome), genome]
vista <- vista[Ntissues<=1]
vista[, label:= fcase(
  heart==1 & limb==0 & midbrain==0, "vista.heart",
  heart==0 & limb==1 & midbrain==0, "vista.limb",
  heart==0 & limb==0 & midbrain==1, "vista.CNS",
  Ntissues==0, "vista.inactive"
)]
dat <- rbind(dat, vista, fill= T)

# Compute kmer
dna <- DNAStringSet(dat$sequence)
kmer_mat <- oligonucleotideFrequency(dna, width = 5)
kmer_freq <- kmer_mat / rowSums(kmer_mat)
kmer_freq <- as.data.table(kmer_freq)
kmer_freq[, label:= dat$label]

# Normalize
agg <- kmer_freq[, lapply(.SD, mean), label]
agg <- as.matrix(agg, 1)
agg <- scale(t(agg))
agg <- agg[apply(agg, 1, function(x) any(x>2)),]
# rownames(agg) <- ifelse(rownames(agg) %in% "GATAA", rownames(agg), "")

# Plot
pdf("pdf/_revision/heatmap_kmer_cluster_validated_enh.pdf", 4)
vl_par(mai= c(.9, .9, .9, 1.5))
vl_heatmap(
  agg,
  cluster.cols = T,
  breaks = seq(-5, 5, length.out= 21),
  legend.title = "kmer freq.\n(z-score)"
)
dev.off()

# pca <- prcomp(t(agg))
pca <- prcomp(scale(kmer_freq[, rownames(agg), with= F]))
pca <- pca$x
rownames(pca) <- kmer_freq$label
labels <- factor(
  rownames(pca),
  c("CNS.validated.enhancer",
    "heart.validated.enhancer",
    "limb.validated.enhancer",
    "vista.CNS",
    "vista.heart",
    "vista.limb",
    "vista.inactive")
)
col <- c("limegreen", "cornflowerblue", "tomato", "limegreen", "cornflowerblue", "tomato", "grey")[labels]
final <- data.table(V1= pca[, "PC1"], V2= pca[, "PC2"], label= labels, col= col)
final[, validated:= grepl("validated", as.character(label))]
setorderv(final, "label", -1)

vl_par()
final[, plot(V1, V2, pch= 21, bg= adjustcolor(col, .5), col= ifelse(validated, "black", NA), xlab= "PC1", ylab= "PC2")]

