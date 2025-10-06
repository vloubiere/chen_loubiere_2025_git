setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
require(orthogene)

# Import clustered TFs ----
dat <- readRDS("Rdata/motif_clusters_paper_3_tissues.rds")
cl <- dat[, .(TF= unlist(tstrsplit(TFs, ","))), cluster]
cl <- na.omit(cl)

# Add all TFs (universe) ----
all_TFs <- readRDS("/groups/stark/vloubiere/motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")$meta
all_TFs <- unique(unlist(all_TFs$TF))
cl <- rbind(cl, data.table(cluster= 0, TF= all_TFs[!all_TFs %in% cl$TF]))

# Find gene IDs ----
mm <- orthogene::map_genes(genes = cl$TF, species = "mouse")
mm <- as.data.table(mm)
mm <- mm[, .SD[1], input]
cl[mm, mouse_id:= i.target, on= "TF==input"]
hs <- orthogene::map_genes(cl[is.na(mouse_id), TF], species = "human")
hs <- as.data.table(hs)
homologs <- orthogene::map_orthologs(hs$target, input_species = "human", output_species = "mouse")
homologs <- as.data.table(homologs)
homologs <- merge(hs[, .(TF= input, human_id= target)], homologs[, .(human_id= input_gene, mouse_id= ortholog_ensg)])
cl[homologs, mouse_id:= i.mouse_id, on= "TF"]
clean <- unique(cl[, .(cluster, mouse_id)])
clean <- na.omit(clean)

# Compute enrichment
genes <- clean[cluster!="0"]
genes[, cluster:= droplevels(cluster)]
genes <- split(genes$mouse_id, genes$cluster)
enr <- vl_GOenrich(geneIDs = genes, geneUniverse.IDs = unlist(unique(clean$mouse_id)), species = "Mm", select = "BP")

# Compute enrichment for housekeeping genes
annot <- readRDS("/groups/stark/nemcko/Pprc1_paper/db/annotations/genes/gene_annotation_Filip_mm10.RDS")
annot <- annot[gene_id %in% clean$mouse_id]
tot.hk <- sum(annot$housekeeping)
hk.enr <- clean[cluster!="0", {
  # Fisher
  hk <-  sum(unique(mouse_id) %in% annot[( housekeeping), gene_id])
  dev <- sum(unique(mouse_id) %in% annot[(!housekeeping), gene_id])
  mat <- c(hk, dev, tot.hk-hk)
  mat <- matrix(c(mat, nrow(annot)-sum(mat)), byrow = T, nrow = 2)
  .(estimate= fisher.test(mat+1)$estimate, pval= fisher.test(mat+1)$p.value)
}, .(cl= cluster)]
hk.enr[, log2OR:= log2(estimate)]
hk.enr[, name:= "housekeeping genes"]

# Select GOs of interest/housekeeping enrichments
sel <- enr[name %in% c("heart development", "limb development", "central nervous system development")]
sel <- rbind(sel, hk.enr, fill= T)
sel[, padj:= p.adjust(pval, "fdr"), cl]
sel[, name:= factor(name, c("heart development", "limb development", "central nervous system development", "housekeeping genes"))]

# Dcast matrices ----
padj <- dcast(sel, name~cl, value.var = "padj")
padj <- as.matrix(padj, 1)
OR <- dcast(sel, name~cl, value.var = "log2OR")
OR <- as.matrix(OR, 1)
# OR[padj>0.05] <- NA
colnames(OR) <- paste0("Cluster ", colnames(OR))
format.padj <- apply(padj, 2, function(x) {
  star <- cut(x, c(-Inf, 1e-5, 1e-3, 1e-2, 5e-2, Inf), c("****\n", "***\n", "**\n", "*\n", "N.S\n"))
  .c <- sapply(formatC(x, format = "e", digits = 1), function(x) paste0("p=", x))
  paste0(star, .c)
})

# Plot ----
vl_par(mai= c(2,2,2,2))
vl_heatmap(
  OR,
  breaks = seq(0, 2.5, 0.1),
  cluster.rows = F,
  show.grid = T,
  show.numbers = format.padj,
  numbers.cex = .4,
  legend.title = "log2OR",
  pdf.file = "pdf/0_paper/GO_enrich_TF_clusters.pdf",
  pdf.cell.size = .4
)
