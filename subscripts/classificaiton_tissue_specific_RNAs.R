setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import annotation ----
if(!exists("annot")) {
  annot <- rtracklayer::import("/groups/stark/shenzhi.chen/projects/accessibility_model_enhancer_design_17112025/annotation/mm10/gencode.vM21.annotation.gtf.gz")
  annot <- as.data.table(annot)
}

# Import RNA data peaks ----
if(!exists("dat")) {
  folder <- "/groups/stark/shenzhi.chen/projects/accessibility_model_enhancer_design_17112025/"
  meta <- readRDS(paste0(folder, "Rdata/mouse_e11.5_ENCODE_20251204/metadata.rds"))
  meta[, files:= paste0(folder, files)]
  meta <- meta[dataset=="RNA"]
  dat <- meta[, {
    fread(files)
  }, .(tissue, rep)]
  # Merge with annotation
  dat <- merge(
    unique(annot[type=="gene", .(gene_id, gene_name, seqnames, start, end, strand)]),
    unique(dat[, .(tissue, rep, gene_id, TPM, FPKM)])
  )
  # Comput emean between replicates
  dat <- dat[, .(log2FPKM= mean(log2(FPKM+0.01))), .(tissue, gene_id, gene_name, seqnames, start, end, strand)]
}

# Look at known tissue-specific genes ----
dat[gene_name=="Nkx2-5"]
dat[gene_name=="Mef2c"]
dat[gene_name=="Gata4"]
dat[gene_name=="Twist1"]
dat[gene_name=="Hoxd9"]
dat[gene_name=="Shh"]
dat[gene_name=="Alb"]
plot(NA, type= "n", xlim= c(-10, 15), ylim= c(0,1))
dat[, lines(density(log2FPKM)), tissue]
abline(v= 0)
setorderv(dat, "log2FPKM", -1)

# Define classes ----
dat[, class:= {
  fcase(
    log2FPKM[1]>0 & log2FPKM[2]<0 & (log2FPKM[1]-log2FPKM[2])>1, paste0(tissue[1], "-specific"),
    log2FPKM[1]>0 & (log2FPKM[1]-log2FPKM[2])>1, paste0(tissue[1], "-enriched"),
    log2FPKM[1]>0, "Shared",
    default= "Inactive"
  )
}, gene_id]

# Clean and save ----
final <- unique(dat[class!="Inactive", .(gene_id, gene_name, seqnames, start, end, strand, class)])
table(final$class)
final[gene_name=="Nkx2-5"]
final[gene_name=="Mef2c"]
final[gene_name=="Gata4"]
final[gene_name=="Twist1"]
final[gene_name=="Hoxd9"]
final[gene_name=="Shh"]
final[gene_name=="Alb"]

saveRDS(final, "Rdata/tissue_specific_genes_RNA.rds")

# # Import K27Ac/K4me1 peaks ----
# mat <- dcast(TPM, gene_id~tissue, value.var = "var")
# mat <- as.matrix(mat, 1)
# # mat <- t(scale(t(mat)))
# mat <- scale(mat)
# mat <- na.omit(mat)
# mat <- mat[apply(mat, 1, function(x) any(x>1)), ]
# mat <- mat[, c("heart", "limb", "forebrain", "midbrain", "hindbrain", "neuralTube", "facialProminence", "liver")]
# cl <- kmeans(mat, 35)$cluster
# clust <- factor(cl)
# agg <- rbindlist(lapply(split(as.data.table(mat), clust), function(x) as.data.table(as.list(apply(x, 2, mean)))))
# agg <- t(scale(t(agg)))
# ord1 <- apply(agg, 1, which.max)
# ord2 <- apply(agg, 1, function(x) any(x>2))
# ord <- as.character(seq(nrow(agg))[order(ord2, ord1, decreasing = T)])
# # ord <- c("6", "4", "23", "20", "24", "25")
# clust <- factor(clust, c(ord, levels(clust)[!levels(clust) %in% ord]))
# 
# vl_par()
# vl_heatmap(
#   mat,
#   cluster.rows= clust,
#   show.row.clusters = "left",
#   breaks= seq(-3, 3, .1)
# )
