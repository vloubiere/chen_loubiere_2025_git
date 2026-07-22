setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

final.pdf <- "pdf/motif_analysis_compar_designed_sequences.pdf"

# Import data ----
folder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/"
# Sequences
sequences <- readRDS("Rdata/final_designed_enhancer_sequences_heart.rds")
sequences[, label:= factor(
  label, 
  c("vista_ts", "rdm", "vista_inactive_all_tissues",
    "evo.ini", "evo.act", "evo.act.acc",
    "ledidi_10_7_ini", "ledidi_10_9_ini", "ledidi_10_12_ini", "ledidi_10_14_ini", "ledidi_12_14_ini",
    "ledidi_10_7", "ledidi_10_9", "ledidi_10_12", "ledidi_10_14", "ledidi_12_14")
)]
# Motifs
mot <- readRDS(paste0(folder, "db/motif/jeff_non-redundant_pwms_jaspar_all.rds"))
# Perc matrices
perc <- readRDS(paste0(folder, "db/motif/jeff_non-redundant_pwms_jaspar_all_percentage.rds"))
# Motif clusters
clust <- readRDS(paste0(folder, "Rdata/final_enhancer_selection/all_motifs_clusters_metatable.rds"))
setnames(clust, c("contrib", "access"), c("enh.contrib", "acc.contrib"))

# Compute motif counts ----
# counts.file <- "db/motifs/motif_counts_selected_sequences_0.00001.rds"
counts.file <- "db/motifs/motif_counts_selected_sequences_0.0001.rds"
if(!file.exists(counts.file)) {
  counts <- vl_motifCounts(sequences = sequences$sequence,
                           pwm_log_odds = mot,
                           genome = "mm10",
                           bg = "genome",
                           p.cutoff = 1e-4)
  # Add info
  counts[, label:= sequences$label]
  counts[, id:= sequences$id]
  counts[, width:= sequences$width]
  counts[, selected:= sequences$selected]
  setcolorder(counts,
              c("label", "id", "width", "selected"))
  # Add Save
  saveRDS(counts,
          counts.file)
}
counts <- readRDS(counts.file)

# Compute motif enrichment based on number of sequences ----
enr.file <- "db/motifs/enrichment_selected_sequences.rds"
if(!file.exists(enr.file)) {
  # Enrichment
  list_counts <- split(counts[, !c("label", "id", "width", "selected")],
                       counts$label)
  enr <- vl_motifEnrich(counts = list_counts[names(list_counts) != "rdm"],
                        control.counts = list_counts[["rdm"]])
  # Add extra columns
  enr[, cl:= factor(cl, levels(sequences$label))]
  enr[, cl:= droplevels(cl)]
  enr[clust, name:= i.cluster, on= "motif"]
  enr[clust, c("enh.contrib", "acc.contrib"):= .(i.enh.contrib, i.acc.contrib), on= "motif"]
  # handle padj of 0
  enr[padj==0, padj:= min(enr[padj>0, padj])]
  # Clean
  enr <- enr[!is.na(acc.contrib) & !is.na(enh.contrib)]
  setcolorder(enr,
              c("cl", "name", "motif", "acc.contrib", "enh.contrib"))
  # Save
  saveRDS(enr, enr.file)
}

# Compute motif enrichment based on motif counts ----
enr2.file <- "db/motifs/counts_enrichment_selected_sequences.rds"
if(!file.exists(enr2.file)) {
  # Melt counts
  enr2 <- data.table::copy(counts)
  enr2 <- melt(enr2,
               id.vars = c("label", "width", "id", "selected"),
               variable.name = "motif")
  # Compute total counts and mean log2
  enr2[, set_hit:= sum(value), .(motif, label)]
  enr2[, mean:= mean(value), .(motif, label)]
  # Normalize for sequence length
  enr2 <- enr2[, .(value= .(value/width)), .(motif, cl= label, set_hit, mean)]
  # Merge
  enr2 <- merge(enr2[cl=="rdm", !"cl"],
                enr2[cl!="rdm"],
                by= "motif")
  enr2[, cl:= droplevels(cl)]
  # Compute statistics
  enr2[, log2_mean_diff:= log2(mean.y/mean.x)]
  enr2[log2_mean_diff > 0, pval:= wilcox.test(value.y[[1]], value.x[[1]], alternative= "greater")$p.value, .(motif, cl)]
  enr2[log2_mean_diff <= 0, pval:= 1]
  # Compute p.adjust
  enr2[, padj:= p.adjust(pval, method = "fdr"), cl]
  # Add extra columns
  enr2[clust, name:= i.cluster, on= "motif"]
  enr2[clust, c("enh.contrib", "acc.contrib"):= .(i.enh.contrib, i.acc.contrib), on= "motif"]
  # Clean
  enr2$set_hit.x <- enr2$mean.x <- enr2$value.x <- enr2$mean.y <- enr2$value.y <- NULL
  setnames(enr2, "set_hit.y", "set_hit")
  setcolorder(enr2,
              c("cl", "name", "motif", "acc.contrib", "enh.contrib"))
  enr2 <- enr2[!is.na(acc.contrib) & !is.na(enh.contrib)]
  # Handle padj of 0
  enr2[padj==0, padj:= min(enr2[padj>0, padj])]
  enr2[, log2OR:= log2_mean_diff]
  # Save
  saveRDS(enr2, enr2.file)
}

# Compute motif enrichment based on number of sequences using only edited counts ----
enr3.file <- "db/motifs/enrichment_selected_sequences_edited.rds"
if(!file.exists(enr.file)) {
  # Split counts and subtract
  list_counts <- split(counts[, !c("label", "id", "width", "selected")],
                       counts$label)
  list_counts$evo.act <- list_counts$evo.act-list_counts$evo.ini
  list_counts$evo.act.acc <- list_counts$evo.act.acc-list_counts$evo.ini
  list_counts$ledidi_10_7 <- list_counts$ledidi_10_7-list_counts$ledidi_10_7_ini
  list_counts$ledidi_10_9 <- list_counts$ledidi_10_9-list_counts$ledidi_10_9_ini
  list_counts$ledidi_10_12 <- list_counts$ledidi_10_12-list_counts$ledidi_10_12_ini
  list_counts$ledidi_10_14 <- list_counts$ledidi_10_14-list_counts$ledidi_10_14_ini
  list_counts$ledidi_12_14 <- list_counts$ledidi_12_14-list_counts$ledidi_12_14_ini
  # Enrichment
  enr3 <- vl_motifEnrich(counts = list_counts[names(list_counts) != "rdm"],
                         control.counts = list_counts[["rdm"]])
  # Add extra columns
  enr3[, cl:= factor(cl, levels(sequences$label))]
  enr3[, cl:= droplevels(cl)]
  enr3[clust, name:= i.cluster, on= "motif"]
  enr3[clust, c("enh.contrib", "acc.contrib"):= .(i.enh.contrib, i.acc.contrib), on= "motif"]
  # handle padj of 0
  enr3[padj==0, padj:= min(enr3[padj>0, padj])]
  # Clean
  enr3 <- enr3[!is.na(acc.contrib) & !is.na(enh.contrib)]
  setcolorder(enr3,
              c("cl", "name", "motif", "acc.contrib", "enh.contrib"))
  # Save
  saveRDS(enr3, enr3.file)
}

# Select best motif per cluster ---
# Fisher method
enr <- readRDS(enr.file)
sel <- enr[cl=="vista_ts", motif[which.min(padj)], name]$V1
enr <- enr[motif %in% sel]
# Motif counts method
enr2 <- readRDS(enr2.file)
sel2 <- enr2[cl=="vista_ts", motif[which.min(padj)], name]$V1
enr2 <- enr2[motif %in% sel2]
# Fisher method edited counts
enr3 <- readRDS(enr3.file)
sel <- enr3[cl=="vista_ts", motif[which.min(padj)], name]$V1
enr3 <- enr3[motif %in% sel]

# Plot top motifs compare parameters ----
top.enr <- 10
col.breaks <- seq(0, 20)
Cc <- c("grey20", "red")

# Plot ----
pdf1 <- tempfile(fileext = ".pdf")
pdf(pdf1,
    width = 12,
    height = 5)
vl_par(mai= c(1,2.5,.9,1.5),
       mfrow= c(1,2))
# Enrichment using fisher test
pl <- plot(enr,
           top.enrich = top.enr,
           cex= .3,
           col= Cc,
           main= "Enrich. vs. rdm genomic seq \nNseq >= 1 motif, fisher",
           color.breaks = col.breaks)

addMotifs(plot.DT = pl,
          pwms = perc,
          cex.width = .5,
          cex.height = .7)
# Enrichment based on motif counts
pl2 <- plot.vl_enr_cl(enr2,
                      top.enrich = top.enr,
                      cex= .35,
                      col= Cc,
                      main= "Enrich. vs. rdm genomic seq \nN motifs / width, wilcox",
                      color.breaks = col.breaks,
                      min.counts = 50)
addMotifs(plot.DT = pl2,
          pwms = perc,
          cex.width = .5,
          cex.height = .7)
# Extract corresponding motif counts ----
sub <- counts[, c("label", "width", unique(pl$motif)), with= FALSE]
sub <- melt(sub,
            id.vars = c("label", "width"),
            variable.name = "motif")
sub[pl, name:= i.name, on= "motif"]
sub[, name:= factor(name, pl[, name[1], keyby= -y]$V1)]
sub <- sub[, .(perc= sum(value>0)/.N,
               mean= mean(value/width)*1001), .(label, name)]
perc.mat <- dcast(sub, name~label, value.var = "perc")
mean.mat <- dcast(sub, name~label, value.var = "mean")
# Balloon plot ----
balloons_plot(size.var = as.matrix(perc.mat, 1),
              color.var = as.matrix(mean.mat, 1),
              col= Cc,
              cex = 2,
              size.legend.title = "% seq >=1 mot.",
              color.legend.title = "Mean mot. counts/kb",
              color.breaks = seq(0, 1.5, length.out= 21))
addMotifs(plot.DT = pl,
          pwms = perc,
          cex.width = .5,
          cex.height = .7)
# Enrichment using fisher test only on edited counts
pl3 <- plot(enr3,
            top.enrich = top.enr,
            cex= .3,
            col= Cc,
            main= "Enrich. vs. rdm genomic seq \nNseq >= 1 motif, fisher",
            color.breaks = col.breaks)
addMotifs(plot.DT = pl3,
          pwms = perc,
          cex.width = .5,
          cex.height = .7)
dev.off()

# Compare enrichment with contrib with different parameters ----
padj.lim <- c(0, 30)
pdf2 <- tempfile(fileext = ".pdf")
pdf(pdf2,
    height = 6.6)
vl_par(mfrow= c(3,3),
       mai= c(.4, .5, .4, .4),
       omi= c(0, 0, 0, .5),
       mgp= c(2, .35, 0))
for(group in c("evo.act", "evo.act.acc", "ledidi_10_7", "ledidi_10_9", "ledidi_10_12", "ledidi_10_14", "ledidi_12_14"))
{
  .c <- enr[cl==group]
  .c[, log10padj:= -log10(padj)]
  Cc <- circlize::colorRamp2(c(padj.lim[1], -log10(0.05), padj.lim[2]),
                             colors = c("grey", "pink", "red"))
  setorderv(.c, "log10padj")
  xlim <- c(-0.001, 0.006)
  ylim <- c(-0.00035, 0.0007)
  .c[, pch:= ifelse(acc.contrib<xlim[2] & enh.contrib<ylim[2], 16, 17)]
  .c[, {
    plot(ifelse(acc.contrib<xlim[2], acc.contrib, xlim[2]),
         ifelse(enh.contrib<ylim[2], enh.contrib, ylim[2]),
         col= Cc(log10padj),
         pch= pch,
         ylab= "Enhancer contrib",
         xlim= xlim,
         ylim= ylim,
         main= group,
         xaxt= "n",
         xlab= NA)
    axis(1, padj= -1.25)
    title(xlab= "Accessibility contrib", line= 1)
    .SD[acc.contrib>0.004 | enh.contrib>0.0004, {
      text(ifelse(acc.contrib<xlim[2], acc.contrib, xlim[2]),
           ifelse(enh.contrib<ylim[2], enh.contrib, ylim[2]),
           name,
           pos= 3,
           cex= .5,
           xpd= T)
    }]
    br <- seq(padj.lim[1], padj.lim[2], length.out= 100)
    par(lwd= .5)
    heatkey(breaks = br,
            col = Cc(br),
            main = "padj (-log10)",
            cex = .6)
    print(".")
  }]
}
dev.off()

# Percentage of enriched motifs with different parameters ----
pdf3 <- tempfile(fileext = ".pdf")
pdf(pdf3,
    width = 10,
    height = 4)
vl_par(mfrow= c(1,2),
       mai= c(1.5, .9, .9, .9),
       cex.main= 1)
pl1 <- enr[, sum(padj<0.05 & log2OR>log2(1.5)), keyby= cl]
vl_barplot(pl1$V1,
           names.arg = pl1$cl,
           ylab= paste0("Number of enriched motifs \ntotal= ", length(unique(enr$name))),
           main= "Fisher method\npadj<0.05, log2OR>log2(1.5)")
pl2 <- enr2[, sum(padj<0.05 & log2OR>log2(1.5)), keyby= cl]
vl_barplot(pl2$V1,
           names.arg = pl2$cl,
           ylab= paste0("Number of enriched motifs \ntotal= ", length(unique(enr2$name))),
           main= "Wilcox method\npadj<0.05, log2OR>log2(1.5)")
dev.off()

# Combine ----
pdftools::pdf_combine(c(pdf1, pdf2, pdf3), final.pdf)
file.show(final.pdf)