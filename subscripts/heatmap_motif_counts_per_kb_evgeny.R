setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Open pdf file ----
pdf.width <- 18
pdf("pdf/heatmap_motif_counts_per_kb_evgeny.pdf",
    width = pdf.width,
    height = 18)
vl_par(mgp= c(2.25, .55, 0.2),
       xaxs= "i",
       lwd= .5,
       cex.lab= .4)
lmat <- imat <- c(1,2,3,4,6,5)
for(i in 1:2)
  lmat <- c(lmat, imat+max(lmat))
layout(matrix(lmat),
       heights = rep(c(.6,.6,.6,.4,.6,5.2), 3))

# For each tissue ----
# for(tiss in "heart"){
for(tiss in c("heart", "limb", "midbrain")){
  # Import sequence info and counts----
  seq.info <- readRDS(paste0("Rdata/final_designed_enhancer_sequences_", tiss, ".rds"))
  counts <- readRDS(paste0("db/motifs/motif_counts_", tiss, "_sequences_0.0001_genome_bg.rds"))
  
  # Select sequences of interest ----
  sel <- c("rdm", "vista_inactive_all_tissues", "vista_ts",
           "ledidi_12_14", "ledidi_12_14_ini")
  lvls <- c("Random genomic seq.", "Inactive VISTA (all tissues)", paste("VISTA", tiss, "enhancers"),
            "LEDIDI design", "Initialization seq.")
  seq.info[, label:= factor(lvls[factor(label, sel)], lvls)]
  counts <- counts[!is.na(seq.info$label),]
  seq.info <- seq.info[!is.na(label)]
  
  # Subset sequences (N vista ts) ----
  Nsubset <- sum(seq.info$label==lvls[3])
  seq.info[, subset:= seq(.N)<=Nsubset, label]
  counts <- counts[seq.info$subset]
  seq.info <- seq.info[(subset), !"subset"]
  
  # Select motif that are enriched in vista_ts/ledidi OR have high contrib ----
  # Import enrichment files
  enr.file <- paste0("db/motifs/motif_enrich_", tiss, "_vs_rdm_seq_fisher.rds")
  enr <- readRDS(enr.file)
  # Selection
  enr[, cl:= factor(lvls[factor(cl, sel)], lvls[3:4])]
  enr <- enr[!is.na(cl) & !is.na(name)]
  if(tiss=="heart") {
    acc.contrib.cutoff <- 0.0015
    enh.contrib.cutoff <- 0.0002
  } else if(tiss=="limb") {
    acc.contrib.cutoff <- 0.0015
    enh.contrib.cutoff <- 0.0002
  }else if(tiss=="midbrain") {
    acc.contrib.cutoff <- 0.001
    enh.contrib.cutoff <- 5e-5
  }
  enr[, enh:= enh.contrib>enh.contrib.cutoff]
  enr[, acc:= acc.contrib>acc.contrib.cutoff]
  enr[, vista:= .SD[cl==lvls[3], padj<0.05 & log2OR>log2(1.5) & set_hit>set_total/20], motif]
  enr[, ledidi:= .SD[cl==lvls[4], padj<0.05 & log2OR>log2(1.5) & set_hit>set_total/20], motif]
  setorderv(enr, c("ledidi", "vista", "enh", "acc", "log2OR"), -1)
  selMot <- enr[(vista | ledidi | enh | acc), .(motif= motif[1]), name]
  top <- enr[selMot$motif, on= "motif"]
  # Order based on ledidi first, then vista_ts enrichments
  top[, ledidi.log2:= log2OR[cl==lvls[4]], motif]
  setorderv(top, c("ledidi", "ledidi.log2"), -1)
  lvl1 <- top[(ledidi), rev(unique(name))]
  top[, vista.log2:= log2OR[cl==lvls[3]], motif]
  setorderv(top, c("vista", "vista.log2"), -1)
  lvl2 <- top[(!ledidi), rev(unique(name))]
  lvls <- rev(unique(rev(c(lvl2, lvl1))))
  top[, name:= factor(name, lvls)]
  setorderv(top, "name")
  top[, motif:= factor(motif, rev(unique(rev(motif))))]
  top.order <- unique(top[, .(name, motif, acc.contrib, enh.contrib)])
  
  # Select counts columns corresponding to enriched motifs ----
  counts <- counts[, unique(top.order$motif), with= FALSE]
  setnames(counts, as.character(top.order$motif), as.character(top.order$name))
  
  # Normalize per kb ----
  counts <- counts/(seq.info$width/1001)
  
  # Retrieve motif enrichment log2OR and padj ----
  log2OR <- dcast(top,
                  cl~name,
                  value.var = "log2OR")
  log2OR <- as.matrix(log2OR, 1)
  padj <- dcast(top,
                cl~name,
                value.var = "padj")
  padj <- as.matrix(padj, 1)
  stars <- padj
  stars[padj<0.05 & log2OR>log2(1.5)] <- "*"
  stars[stars!="*"] <- NA
  
  # Title ----
  par(mai= c(0,0,0,0))
  plot.new()
  text(.5, .5, tiss, cex= 2)
  
  # Initialize plotting ----
  adj <- (pdf.width-(0.1*ncol(counts)))/2
  par(mai= c(.05, adj, .05, adj))
  
  # Plot contributions (barplots) ----
  vl_barplot(top.order$enh.contrib,
             xaxt= "n",
             ylab= "Enh.contrib")
  vl_barplot(top.order$acc.contrib,
             xaxt= "n",
             ylab= "Acc.contrib")
  
  # Heatmap enrichments ----
  vl_heatmap(log2OR,
             cluster.rows = F,
             show.numbers = stars,
             legend.title = "log2OR",
             breaks = seq(-3, 3, length.out= 21),
             legend.cex = .8)
  
  # Heatmap parameters ----
  br <- seq(0, 3, length.out= 21)
  Cc <- colorRampPalette(c("white", "red"))(length(br))
  leg <- "Motif/kb"
  
  # Plot normalized motif counts ----
  par(mai= c(.5, adj, .05, adj))
  vl_heatmap(counts,
             cluster.rows = seq.info$label,
             show.rownames = F,
             breaks = br,
             col= Cc,
             legend.title = leg,
             show.row.clusters = "left",
             legend.cex = .8)
  par(mai= c(0,0,0,0))
  plot.new()
}
dev.off()