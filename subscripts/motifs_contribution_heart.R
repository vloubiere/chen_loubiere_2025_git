setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
library(ggplot2)
library(ggrepel)

# Import data ----
folder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/"
clust <- readRDS(paste0(folder, "Rdata/final_enhancer_selection/all_motifs_clusters_metatable.rds"))
setnames(clust, c("contrib", "access"), c("enh.contrib", "acc.contrib"))
dat <- clust[, c(list(motif= motif[which.max(enh.contrib)]),
                 lapply(.SD, max)), cluster, .SDcols= c("acc.contrib", "enh.contrib")]
enr <- readRDS("db/motifs/enrichment_selected_sequences.rds") 
dat[enr[cl=="vista_ts"], `-log10(padj) motif enrich.`:= -log10(i.padj), on= "motif"]

# Filtrer les points à annoter
to_label <- dat[acc.contrib > 0.002 | enh.contrib > 0.0004]

# Plot minimaliste mais lisible
ggplot(dat, aes(x = acc.contrib, y = enh.contrib)) +
  geom_point(size = 1, color = "black") +  # Tous les points en noir
  geom_point(
    data = to_label,
    aes(x = acc.contrib, y = enh.contrib),
    color = "red",
    size = 1.5
  ) +  # Points à annoter en rouge
  geom_text_repel(
    data = to_label,
    aes(label = cluster),
    size = 3,
    color = "red",  # Noms en rouge
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.2
  ) +
  labs(x = "Accessibility contribution", y = "Enhancer contribution") +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.background = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05))

# Plot: -log10(padj) motif enrich. vs enhancer contribution
to_label2 <- dat[enh.contrib > 0.0004 | `-log10(padj) motif enrich.` > 7.5]
ggplot(dat, aes(x = enh.contrib, y = `-log10(padj) motif enrich.`)) +
  geom_point(size = 1, color = "black") +
  geom_point(
    data = to_label2,
    aes(x = enh.contrib, y = `-log10(padj) motif enrich.`),
    color = "red",
    size = 1.5
  ) +
  geom_text_repel(
    data = to_label2,
    aes(label = cluster),
    size = 3,
    color = "red",
    max.overlaps = Inf,
    box.padding = 0.3,
    point.padding = 0.2,
    segment.size = 0.2
  ) +
  labs(
    x = "Enhancer contribution",
    y = expression(-log[10]*"(padj) motif enrichment")
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_text(),
    axis.text = element_text(),
    plot.background = element_blank(),
    panel.border = element_blank()
  ) +
  scale_x_continuous(expand = expansion(mult = 0.05)) +
  scale_y_continuous(expand = expansion(mult = 0.05))