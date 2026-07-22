setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Open pdf file ----
pdf.width <- 18
pdf("pdf/heatmap_motif_counts_per_kb.pdf",
    width = pdf.width,
    height = 6.9)

# For each tissue ----
# for(tiss in "heart"){
for(tiss in c("heart", "limb", "midbrain")){
  # Import sequence info and counts----
  seq.info <- readRDS(paste0("Rdata/final_designed_enhancer_sequences_", tiss, ".rds"))
  counts <- readRDS(paste0("db/motifs/motif_counts_", tiss, "_sequences_0.0001_genome_bg.rds"))
  
  # Select sequences of interest ----
  sel <- c("rdm", "vista_inactive_all_tissues", "vista_ts",
           "ledidi_12_14_ini", "evo.act", "evo.act.acc", "ledidi_12_14")
  seq.info[, label:= factor(label, sel)]
  counts <- counts[!is.na(seq.info$label),]
  seq.info <- seq.info[!is.na(label)]
  
  # Subset sequences for first heatmap ----
  Nsubset <- sum(seq.info$label=="vista_ts")
  seq.info[, subset:= seq(.N)<=Nsubset, label]
  
  # Select motif that are enriched in vista_ts/ledidi OR have high contrib ----
  # Import enrichment files
  enr.file <- paste0("db/motifs/motif_enrich_", tiss, "_vs_rdm_seq_fisher.rds")
  enr <- readRDS(enr.file)
  enr[, cl:= factor(cl, c("vista_ts", "ledidi_12_14"))]
  enr <- enr[!is.na(cl) & !is.na(name)]
  # Selection
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
  enr[, vista:= .SD[cl=="vista_ts", padj<0.05 & log2OR>log2(1.5) & set_hit>set_total/20], motif]
  enr[, ledidi:= .SD[cl=="ledidi_12_14", padj<0.05 & log2OR>log2(1.5) & set_hit>set_total/20], motif]
  setorderv(enr, c("ledidi", "vista", "enh", "acc", "log2OR"), -1)
  selMot <- enr[(vista | ledidi | enh | acc), .(motif= motif[1]), name]
  top <- enr[selMot$motif, on= "motif"]
  # Order based on ledidi first, then vista_ts enrichments
  top[, ledidi.log2:= log2OR[cl=="ledidi_12_14"], motif]
  setorderv(top, c("ledidi", "ledidi.log2"), -1)
  lvl1 <- top[(ledidi), rev(unique(name))]
  top[, vista.log2:= log2OR[cl=="vista_ts"], motif]
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
  
  # Subset vista and ledidi sequences ----
  sub.info <- seq.info[label %in% c("vista_ts", "ledidi_12_14")]
  sub.counts <- counts[seq.info$label %in% c("vista_ts", "ledidi_12_14")]
  sub.info[, annot:= fcase(
    label=="vista_ts", "vista_ts",
    selected, "Ledidi: passed filtering",
    default = "Ledidi: did not pass filtering"
  )]
  
  # Kmeans clustering ----
  scaled <- t(scale(t(log2(sub.counts+.5))))
  scaled[is.na(scaled)] <- 0
  scaled <- as.data.table(scaled)
  rownames(scaled) <- sub.info$id
  set.seed(1)
  sub.info$cl <- kmeans(scaled, centers = 8)$cluster
  sub.info[, main.group:= names(rev(sort(table(annot))))[1], cl]
  sub.info[, cl.name:= paste0("Cluster ", cl)]
  sub.info[, line.idx:= .I]
  
  # For each cluster, find sequence with closest distance ----
  NselPerClust <- 2
  sub.info.sel <- sub.info[(selected)]
  sub.counts.sel <- sub.counts[sub.info$selected,]
  sub.info.sel[, best:= {
    .c <- scaled[line.idx,]
    .center <- matrix(apply(.c, 2, mean), nrow= 1)
    colnames(.center) <- colnames(.c)
    .d <- stats::dist(rbind(.center, .c))
    .d <- as.matrix(.d)
    best <- id[order(.d[-1,1])]
    if(length(best)>NselPerClust)
      best <- best[1:NselPerClust]
    id %in% best
  }, cl]
  
  # Mean all seq per cluster ---
  m.counts.all <- as.data.table(sub.counts)[, lapply(.SD, mean), keyby= .(cl= sub.info$cl)]
  m.counts.all <- as.matrix(m.counts.all, 1)
  rownames(m.counts.all) <- paste0("Cl.", rownames(m.counts.all), " (all seq.)")
  annot1 <- sub.info[, main.group[1], keyby= cl]$V1
  
  # Mean seq that passed filtering per cluster ----
  m.counts.sel <- as.data.table(sub.counts.sel)[, lapply(.SD, mean), keyby= .(cl= sub.info.sel$cl)]
  m.counts.sel <- as.matrix(m.counts.sel, 1)
  rownames(m.counts.sel) <- paste0("Cl.", rownames(m.counts.sel), " (passed filtering)")
  annot2 <- sub.info.sel[, main.group[1], keyby= cl]$V1
  
  # Select relevant motifs ----
  occuringMot <- apply(m.counts.sel, 2, function(x) max(x)>.5)
  
  # Combine the two ----
  cmb <- rbind(m.counts.all, m.counts.sel)
  cmb <- cmb[, (occuringMot)]
  cmb.scaled <- t(scale(t(log2(cmb+.5))))
  annot <- c(annot1, annot2)
  
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
  
  # Plot different heatmaps ----
  for(hm in c("groups", "clustering", "mean.counts")) {
    if(hm!="mean.counts") {
      # Initialize plotting ----
      adj <- (pdf.width-(0.1*ncol(counts)))/2
      # vl_par(mai= c(.05, 1.25, .05, 1),
      vl_par(mai= c(.05, adj, .05, adj),
             omi= c(.5, 0, .5, 0),
             mgp= c(2.25, .55, 0.2),
             xaxs= "i",
             lwd= .5,
             cex.lab= .45)
      layout(matrix(c(1,2,3,5,4)),
             heights = c(.4,.4,.4,.6,6))
    } else if(hm=="mean.counts") {
      # Initialize plotting ----
      adj <- (pdf.width-(0.1*sum(occuringMot)))/2
      # par(mai= c(.05, 1.5, .05, 2))
      par(mai= c(.05, adj, .05, adj),
          cex.lab= .45)
      layout(matrix(c(1,2,3,7,4,5,6)),
             heights = c(.4,.4,.4,.6,2,2,2))
      # Subset motifs ----
      top.order <- top.order[occuringMot]
      log2OR <- log2OR[, occuringMot]
      stars <- stars[, occuringMot]
    }
    
    # Plot contributions (barplots) ----
    vl_barplot(top.order$enh.contrib,
               xaxt= "n",
               ylab= "Enh.contrib")
    vl_barplot(top.order$acc.contrib,
               xaxt= "n",
               ylab= "Acc.contrib")
    mtext(text = tiss, outer = T, cex= 2)
    
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
    if(hm=="groups") {
      vl_heatmap(counts[seq.info$subset,],
                 cluster.rows = seq.info[(subset),label],
                 show.rownames = F,
                 breaks = br,
                 col= Cc,
                 legend.title = leg,
                 show.row.clusters = "left",
                 legend.cex = .8)
      plot.new()
    } else if(hm=="clustering") {
      vl_heatmap(sub.counts,
                 cluster.rows = sub.info$cl.name,
                 row.annotations = sub.info$annot,
                 show.rownames = F,
                 breaks = br,
                 col= Cc,
                 legend.title = leg,
                 row.annotations.title = "Sequence category",
                 show.row.clusters = "left",
                 legend.cex = .8)
      plot.new()
    } else if(hm=="mean.counts") {
      # Plot mean counts
      vl_heatmap(cmb,
                 cluster.rows = cmb.scaled,
                 row.annotations = annot,
                 show.rownames = T,
                 show.colnames = F,
                 col= Cc,
                 breaks = br,
                 legend.title = "Mean mot. counts/kb",
                 row.annotations.title = "Most represented seq. category in cluster",
                 legend.cex = .8,
                 show.numbers = round(cmb, 1),
                 numbers.cex = .45)
      # Final selection 
      m.counts.final <- as.data.table(sub.counts.sel)[(sub.info.sel$best)]
      m.counts.final <- as.matrix(m.counts.final)
      m.counts.final <- m.counts.final[, (occuringMot)]
      rownames(m.counts.final) <- sub.info.sel[(best), paste0("Cluster ", cl, " (id=", id, ")")]
      annot3 <- sub.info.sel[(best), rep(main.group[1], .N), keyby= cl]$V1
      cl <- vl_heatmap(m.counts.final,
                       cluster.rows = t(scale(t(log2(m.counts.final+.5)))),
                       row.annotations = annot3,
                       show.rownames = T,
                       col= Cc,
                       breaks = br,
                       legend.title = "Mean mot. counts/kb",
                       row.annotations.title = "Most represented seq. category in cluster",
                       legend.cex = .8,
                       show.numbers = m.counts.final,
                       numbers.cex = .6)
      cl$rows[, id:= tstrsplit(name, "id=|)", keep= 2)]
      cl$rows[, annot:= annot3]
      plot.new()
      plot.new()
      
      # Add best sequence to object ----
      seq.info.file <- paste0("Rdata/final_designed_enhancer_sequences_", tiss, ".rds")
      saveBest <- readRDS(seq.info.file)
      saveBest$selected_validation <- NULL
      saveBest[cl$rows, selected_validation:= i.y.pos, on= "id"]
      saveBest[cl$rows, selected_validation_name:= paste0(i.name, " Cluster majority -> ", i.annot), on= "id"]
      saveRDS(saveBest, seq.info.file)
      
      # Save motifs ----
      cl$cols[top, motif:= i.motif, on= "name"]
      saveRDS(cl$cols, paste0("Rdata/final_designed_enhancer_motifs_", tiss, ".rds"))
    }
  }
}
dev.off()