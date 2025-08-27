setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motif clusters ----
mot <- readRDS("Rdata/annotated_PWMs.rds")

# Open pdf file ----
pdf.width <- 15
pdf("pdf/heatmap_clustering_filtered_sequences_evgeny.pdf",
    width = pdf.width,
    height = 3*3)
vl_par(lwd= .5,
       cex.lab= .4,
       mfrow= c(3, 1))

# For each tissue ----
# for(tiss in "heart"){
for(tiss in c("heart", "limb", "midbrain")){
  # Import sequence info and counts----
  seq.info <- readRDS(paste0("Rdata/final_designed_enhancer_sequences_", tiss, ".rds"))
  counts <- readRDS(paste0("db/motifs/motif_counts_", tiss, "_sequences_0.0001_subject_bg.rds"))
  
  # Subtract initiation sequences ----
  counts <- counts[seq.info$label=="ledidi_12_14"]-counts[seq.info$label=="ledidi_12_14_ini"]
  seq.info <- seq.info[label=="ledidi_12_14"]
  
  # Keep only selected sequences ----
  counts <- counts[seq.info$selected]
  seq.info <- seq.info[(selected)]
  
  # Select motifs that are implanted in at least 50 sequences ----
  counts[, id:= seq.info$id]
  .m <- melt(counts, id.vars = "id", variable.name = "motif")
  .m[mot, cluster:= i.cluster, on= "motif"]
  .m[, Npos:= sum(value>0), motif]
  .m[, mean:= mean(value), motif]
  sel <- .m[!is.na(cluster) & Npos>50 & mean>0.1, .(motif= motif[which.max(mean)]), cluster]
  counts <- dcast(.m[motif %in% sel$motif], id~cluster, value.var = "value")
  counts <- as.matrix(counts, 1)
  counts <- counts[, apply(counts, 2, function(x) any(x>0))]
  
  # Clustering based on counts (no normalization) ----
  km <- kmeans(counts, centers = 12)
  sel.seq <- data.table(id= character(),
                        name= character())
  for(i in seq(nrow(km$centers))) {
    idx <- which(km$cluster==i)
    .d <- dist(
      rbind(km$centers[i,],
            counts[idx, ])
    )
    .d <- as.matrix(.d)
    .d <- .d[-1,1]
    closest.seq.id <- names(.d)[which.min(.d)]
    sel.seq <- rbind(
      sel.seq,
      data.table(
        id= closest.seq.id,
        name= paste0(closest.seq.id, " (cluster ", i, "; n= ", length(idx), ")")
      )
    )
  }
  sel.seq[, label:= "ledidi_12_14"]
  
  # Subset counts ----
  sub.counts <- counts[sel.seq$id,]
  rownames(sub.counts) <- sel.seq$name
  sub.counts <- sub.counts[, order(apply(sub.counts, 2, function(x) sum(x>0)))]
  sub.counts <- sub.counts[, apply(counts, 2, function(x) any(x>0))]
  
  # Initialize plotting ----
  adj <- (pdf.width-(0.1*ncol(sub.counts)))/2
  par(mai= c(.9, adj, .9, adj))
  
  # Heatmap ----
  br <- seq(-1, 4, .1)
  Cc <- c("cornflowerblue", "white", "tomato")
  col <- circlize::colorRamp2(c(min(br), 0, max(br)), colors = Cc)(br)
  cl <- vl_heatmap(sub.counts,
                   breaks = br,
                   col= col,
                   legend.title = "Motif counts (designed-initialization)",
                   legend.cex = .8,
                   show.numbers = T,
                   numbers.cex= .5,
                   main= paste(tiss, "selected sequences (closest distance to kmeans centroids, k= 12)"))
  
  # Add selection to seq info ---
  sel.seq[, name:= factor(name, cl$rows[order(y.pos), name])]
  seq.info.file <- paste0("Rdata/final_designed_enhancer_sequences_", tiss, ".rds")
  seq.info.update <- readRDS(seq.info.file)
  seq.info.update$selected_validation_name <- NULL
  seq.info.update[sel.seq, selected_validation_name:= i.name, on= c("label", "id")]
  saveRDS(seq.info.update, seq.info.file)
  
  # Save selected motifs ----
  sel <- sel[rev(colnames(sub.counts)), on= "cluster"]
  saveRDS(sel, 
          paste0("Rdata/final_designed_enhancer_motifs_", tiss, ".rds"))
}
dev.off()