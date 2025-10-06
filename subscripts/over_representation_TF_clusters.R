setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
require(metap)
require(AUCell)
require(Matrix)

# Import data ----
dat <- readRDS("db/single_cell/subsetted_sc_dataset.rds")

# Import motifs ----
mot <- readRDS("Rdata/motif_clusters_paper_3_tissues_single_TF.rds")
mot <- mot[, .(TF= unlist(TFs)), .(name, cluster)]
TFs <- unique(unlist(mot$TF))
TFs[TFs=="HXA13"] <- "HOXA13"
TFs[TFs=="HXC9"] <- "HOXC9"
# TFs[TFs=="POU5F1B"] <- "POU5F1"
TFs <- intersect(TFs, rownames(dat$counts))
sub <- dat$counts[TFs,]

# Compute enrichment ----
enr.file <- "db/motifs/fisher_test_TF_per_cell_type_sc.rds"
if(!file.exists(enr.file)) {
  enr <- lapply(seq(nrow(sub)), function(i) {
    dat$cells[, {
      fisher.test(sub[i,] > 0, dat$cells$cluster==cluster)[c("estimate", "p.value")]
    }, cluster]
  })
  names(enr) <- rownames(sub)
  res <- rbindlist(enr, idcol= "gene_name")
  res[, padj:= p.adjust(p.value, "fdr"), cluster]
  res[, log2OR:= ifelse(estimate==0, log2(min(estimate[estimate>0])), log2(estimate))]
  saveRDS(res, enr.file)
} else
  res <- readRDS(enr.file)

res[mot, TF.cluster:= i.cluster, on= "gene_name==TF"]

test <- res[cluster!="Other" & TF.cluster!="0"]
test <- test[gene_name %in% res[log2OR>1 & padj<1e-5, gene_name]]
test[, TF.cluster:= factor(TF.cluster, c("4", "1", "3", "2"))]
test[, cluster:= factor(cluster, c("Brain", "Neural tube", "Cardiac muscle lineages", "Limb mesenchyme"))]
mat <- dcast(test, TF.cluster+gene_name~cluster, value.var = "log2OR")

vl_par()
vl_heatmap(as.matrix(mat[, -1], 1), cluster.rows = mat$TF.cluster)
