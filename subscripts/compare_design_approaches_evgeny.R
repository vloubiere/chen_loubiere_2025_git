setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
mot <- readRDS("Rdata/annotated_PWMs.rds")

# Open pdf file ----
pdf("pdf/compare_motif_enrichment_different_design_approaches_evgeny.pdf",
    width = 9.5,
    height = 20)
vl_par(mfcol= c(4, 3))

# For each tissue ----
for(tiss in c("heart", "limb", "midbrain")){
  # Selected labels ----
  sel <- c("vista_ts", "vista_inactive_all_tissues",
           "ledidi_12_14_ini", "evo.act.acc", "ledidi_12_14")
  lvls <- c(paste("VISTA", tiss, "enhancers"), "Inactive VISTA (all tissues)",
            "Initialization seq.", "EVO design (act+acc)", "LEDIDI design")
  
  # Import enrichment files ----
  enr.file <- paste0("db/motifs/motif_enrich_", tiss, "_vs_rdm_seq_fisher.rds")
  enr <- readRDS(enr.file)
  enr[, cl:= factor(lvls[factor(cl, sel)], lvls)]
  enr <- enr[!is.na(cl) & !is.na(name)]
  
  # Select top motifs ----
  top <- enr[, {
    # If any enrichment in vista ts
    if(any(cl==lvls[1] & padj<0.05 & log2OR>log2(1.5))) {
      .SD[cl==lvls[1] & log2OR>log2(1.5)][which.min(padj), motif]
      
    # If any enrichment 
    } else if (any(padj<0.05 & log2OR>log2(1.5))) {
      .SD[log2OR>log2(1.5)][which.min(padj), motif]
      
    # Else
    } else
      NA_character_
  }, name]$V1
  top <- enr[motif %in% na.omit(top)]
  
  # Order and cap ----
  top[, padj:= padj+1e-100]
  
  # Plot motif enrichment ----
  par(mai= c(1, 1.8, .9, .8))
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
  par(mai= c(1.8, .9, 1.8, .9),
      mgp= c(2, .35, 0))
  xlim <- c(-0.001, 0.005)
  ylim <- c(-0.0002, 0.0006)
  setorderv(top, "padj", -1)
  top[cl %in% c("EVO design (act+acc)", "LEDIDI design"), {
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
           cex= .4,
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