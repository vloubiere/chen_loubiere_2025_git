setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

final.pdf <- "pdf/motif_analysis_best_ledidi_parameters.pdf"

# Select best param ----
sel.param <- c("vista_ts",
               "vista_inactive_all_tissues",
               "evo.ini",
               "evo.act",
               "evo.act.acc",
               "ledidi_12_14_ini",
               "ledidi_12_14")
best.ledidi <- "ledidi_12_14"
# Import perc matrices ----
s.folder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/"
perc <- readRDS(paste0(s.folder, "db/motif/jeff_non-redundant_pwms_jaspar_all_percentage.rds"))

# Import motif counts ----
counts <- readRDS("db/motifs/motif_counts_selected_sequences_0.0001.rds")
counts <- counts[label %in% sel.param]
counts[, label:= droplevels(label)]

# Import enrichments and select best motif per cluster ---
# Fisher method
enr <- readRDS("db/motifs/enrichment_selected_sequences.rds")
sel <- enr[cl=="vista_ts", motif[which.min(padj)], name]$V1
enr <- enr[motif %in% sel & cl %in% sel.param]
enr[, cl:= droplevels(cl)]
# Motif counts method
enr2 <- readRDS("db/motifs/counts_enrichment_selected_sequences.rds")
sel2 <- enr2[cl=="vista_ts", motif[which.min(padj)], name]$V1
enr2 <- enr2[motif %in% sel2 & cl %in% sel.param]
enr2[, cl:= droplevels(cl)]
# Motif counts method
enr3 <- readRDS("db/motifs/enrichment_selected_sequences_edited.rds")
sel3 <- enr3[cl=="vista_ts", motif[which.min(padj)], name]$V1
enr3 <- enr3[motif %in% sel3 & cl %in% sel.param]
enr3[, cl:= droplevels(cl)]

# Plot top motifs compare parameters ----
top.enr <- 10
col.breaks <- seq(0, 20)
Cc <- c("grey20", "red")

# Plot ----
pdf1 <- tempfile(fileext = ".pdf")
pdf(pdf1,
    width = 10,
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
                      cex= .5,
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
    height = 2.2)
vl_par(mfrow= c(1,3),
       mai= c(.4, .5, .4, .4),
       omi= c(0, 0, 0, .5),
       mgp= c(2, .35, 0))
for(group in c("evo.act", "evo.act.acc", best.ledidi))
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
    width = 6,
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