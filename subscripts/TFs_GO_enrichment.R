setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
require(orthogene)

# Import clustered TFs ----
TFs <- readRDS("Rdata/motif_clusters_paper_3_tissues.rds")
TFs <- TFs[, .(mouse_id= unlist(mouse_id)), cluster]
TFs <- unique(na.omit(TFs))
universe <- readRDS("/groups/stark/vloubiere/motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")
universe <- na.omit(unlist(universe$meta$mouse_id))

# Compute GO enrichment ----
enr <- vl_GOenrich(geneIDs = split(TFs$mouse_id, TFs$cluster),
                   geneUniverse.IDs = universe,
                   species = "Mm",
                   select = "BP")

# Compute enrichment for housekeeping genes ----
annot <- readRDS("/groups/stark/nemcko/Pprc1_paper/db/annotations/genes/gene_annotation_Filip_mm10.RDS")
annot <- annot[gene_id %in% universe, .(gene_id, housekeeping)]
annot <- unique(annot)
tot <- nrow(annot)
tot.hk <- sum(annot$housekeeping)
hk.enr <- TFs[, {
  # Count
  hk.pos <- sum( mouse_id %in% annot$gene_id[annot$housekeeping])
  hk.neg <- sum(!mouse_id %in% annot$gene_id[annot$housekeeping])
  # Contingency
  mat <- c(hk.pos, hk.neg, tot.hk-hk.pos)
  mat <- c(mat, tot-sum(mat))
  mat <- matrix(mat, byrow = T, nrow = 2)
  # Fisher
  .(
    estimate= fisher.test(mat+1, alternative = "greater")$estimate,
    pval= fisher.test(mat+1, alternative = "greater")$p.value
  )
}, .(cl= cluster)]
hk.enr[, log2OR:= log2(estimate)]
hk.enr[, name:= "housekeeping genes"]

# Select GOs of interest/housekeeping enrichments ----
sel <- enr[name %in% c("heart development", "limb development", "central nervous system development")]
sel <- rbind(sel, hk.enr, fill= T)
sel[, padj:= p.adjust(pval, "fdr"), cl]
sel[, name:= factor(name, c("heart development", "limb development", "central nervous system development", "housekeeping genes"))]

# Dcast matrices ----
padj <- dcast(sel, cl~name, value.var = "padj")
padj <- as.matrix(padj, 1)
OR <- dcast(sel, cl~name, value.var = "log2OR")
OR <- as.matrix(OR, 1)
OR[padj>0.05] <- NA
rownames(OR) <- paste0("Cluster ", rownames(OR))
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
