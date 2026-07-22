setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
mot <- readRDS("Rdata/annotated_PWMs.rds")

# Open pdf file ----
pdf("pdf/compare_motif_enrichment_different_design_approaches.pdf",
    width = 12,
    height = 20)
vl_par(mfcol= c(4, 3))

# For each tissue ----
for(tiss in c("heart", "limb", "midbrain")){
  # Selected labels ----
  sel <- c("vista_ts", "rdm", "vista_inactive_all_tissues",
           "evo.ini", "evo.act", "evo.act.acc",
           "ledidi_10_7_ini", "ledidi_10_9_ini", "ledidi_10_12_ini", "ledidi_10_14_ini", "ledidi_12_14_ini",
           "ledidi_10_7", "ledidi_10_9", "ledidi_10_12", "ledidi_10_14", "ledidi_12_14")
  
  # Import enrichment files ----
  enr.file <- paste0("db/motifs/motif_enrich_", tiss, "_vs_rdm_seq_fisher.rds")
  enr <- readRDS(enr.file)
  enr[, cl:= factor(cl, sel)]
  enr <- enr[!is.na(cl)]
  
  # Select top motifs ----
  top <- enr[cl=="vista_ts" & !is.na(name) & log2OR>log2(1.5), motif[which.min(padj)], name]$V1
  top <- enr[motif %in% top]
  
  # Order and cap ----
  top[, uncapped.log2OR:= log2OR]
  top[, log2OR:= ifelse(log2OR>3, 3, log2OR)]
  top[padj==0, padj:= min(top$padj[top$padj>0])]
  
  # Plot motif enrichment ----
  par(mai= c(1, 1.75, .9, .7))
  pl <- plot(top,
             top.enrich = 12,
             cex= .8,
             col= c("salmon2", "red"),
             main= tiss,
             color.breaks = seq(0, 10))
  addMotifs(plot.DT = pl,
            pwms = mot$pwm.perc,
            cex.width = .5,
            cex.height = .7)
  
  # Scatterplot contrib & enrichment ----
  par(mai= c(1.6, 1.1, 1.6, 1.1),
      mgp= c(2, .35, 0))
  xlim <- c(-0.001, 0.005)
  ylim <- c(-0.0002, 0.0006)
  setorderv(top, "padj", -1)
  top[cl %in% c("evo.act.acc", "ledidi_12_14"), {
    plot(ifelse(acc.contrib<xlim[2], acc.contrib, xlim[2]),
         ifelse(enh.contrib<ylim[2], enh.contrib, ylim[2]),
         col= ifelse(padj<0.05, "red", "grey"),
         pch= ifelse(acc.contrib<xlim[2] & enh.contrib<ylim[2], 16, 17),
         ylab= "Enhancer contrib",
         xlim= xlim,
         ylim= ylim,
         xaxt= "n",
         xlab= NA,
         main= cl)
    axis(1, padj= -1.25)
    title(xlab= "Accessibility contrib", line= 1)
    # Add names
    .SD[acc.contrib>0.002 | enh.contrib>0.0002 | padj<0.05, {
      text(ifelse(acc.contrib<xlim[2], acc.contrib, xlim[2]),
           ifelse(enh.contrib<ylim[2], enh.contrib, ylim[2]),
           name,
           col= ifelse(padj<0.05, "red", "black"),
           pos= 3,
           cex= .5,
           xpd= T)
    }]
  }, cl]
  
  # Barplot number of enriched motifs ----
  pl1 <- enr[, length(unique(name[padj<0.05 & log2OR>log2(1.5)])), keyby= cl]
  par(mai= c(2, 1.1, 2, 1.1))
  vl_barplot(pl1$V1,
             names.arg = pl1$cl,
             ylab= paste0("Number of enriched motifs \ntotal= ", length(unique(enr$name))),
             main= "Fisher method\npadj<0.05, log2OR>log2(1.5)")
}
dev.off()