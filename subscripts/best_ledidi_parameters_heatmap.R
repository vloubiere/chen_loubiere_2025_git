setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Select best param ----
sel.param <- c("rdm",
               "vista_inactive_all_tissues",
               "vista_ts",
               "ledidi_12_14_ini",
               "evo.act",
               "evo.act.acc",
               "ledidi_12_14")
best.ledidi <- "ledidi_12_14"

# Import enrichments and select clusters ---
enr <- readRDS("db/motifs/enrichment_selected_sequences.rds")
enr <- enr[cl %in% c("vista_ts", best.ledidi)]
enr[, cl:= droplevels(cl)]
contrib <- unique(enr[, .(name, acc.contrib, enh.contrib)])

# Import prediction scores heart enhancer test dataset ----
pred <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/final_enhancer_selection/vista_test_dataset_predicted_acc_act_heart.rds")

# Set predicted act and contrib cutoffs ----
predict.act.cutoff <- .8
acc.contrib.cutoff <- 0.0015
enh.contrib.cutoff <- 0.0002
pdf("pdf/cutoffs_selected_enhancers.pdf",
    width = 1.5,
    height = 4.5)
vl_par(mfrow= c(3,1),
       mai= c(.2,.7,.2,.2),
       mgp= c(2.5, .5, 0),
       omi= c(0, 0, .4, .2))
# Predicted act
vl_boxplot(pred$vista_act,
           violin= T,
           outline= T,
           ylab= "Predicted act.",
           main= "Vista ts enh.")
abline(h= predict.act.cutoff)
# Contrib cutoff
vl_boxplot(contrib$acc.contrib,
           violin= T,
           ylab= "Mean acc. contrib/motif",
           main= "All motifs")
abline(h= acc.contrib.cutoff)
# Enhancer cutoff
vl_boxplot(contrib$enh.contrib,
           violin= T,
           ylab= "Mean enh. contrib/motif",
           main= "All motifs")
abline(h= enh.contrib.cutoff)
dev.off()

# Select motifs that are enriched or have high contribution scores ----
enr[, vista:= .SD[cl=="vista_ts", padj<0.05 & log2OR>log2(1.5) & set_hit>set_total/20], motif]
enr[, ledidi:= .SD[cl==best.ledidi, padj<0.05 & log2OR>log2(1.5) & set_hit>set_total/20], motif]
enr[, enh:= enh.contrib>enh.contrib.cutoff]
enr[, acc:= acc.contrib>acc.contrib.cutoff]
setorderv(enr, c("ledidi", "vista", "enh", "acc", "log2OR"), -1)
selMot <- enr[cl %in% c("vista_ts", best.ledidi) & (vista | ledidi | enh | acc), .(motif= motif[1]), name]
enr <- enr[selMot$motif, on= "motif"]

# Order based on ledidi first, then vista_ts enrichments ----
enr[, ledidi.log2:= log2OR[cl==best.ledidi], motif]
setorderv(enr, c("ledidi", "ledidi.log2"), -1)
lvl1 <- enr[(ledidi), rev(unique(name))]
enr[, vista.log2:= log2OR[cl=="vista_ts"], motif]
setorderv(enr, c("vista", "vista.log2"), -1)
lvl2 <- enr[(!ledidi), rev(unique(name))]
lvls <- rev(unique(rev(c(lvl2, lvl1))))
enr[, name:= factor(name, lvls)]
setorderv(enr, "name")
enr[, motif:= factor(motif, rev(unique(rev(motif))))]

# Import motif counts and select relevant sequences----
counts <- readRDS("db/motifs/motif_counts_selected_sequences_0.0001.rds")
counts <- counts[label %in% sel.param, c("label", "id", "width", "selected", levels(enr$motif)), with= FALSE]
counts[, label:= factor(label, sel.param)]
setnames(counts, levels(enr$motif), levels(enr$name))
# Normalize per kb
cols <- levels(enr$name)
counts[, (cols):= lapply(.SD, function(x) x/width*1001), width, .SDcols= cols]
# Split Seq info and counts
seq.info <- counts[, .(label, id, width, selected)]
counts <- counts[, levels(enr$name), with= FALSE]

# Remove motifs with very high counts (~1/50bp)
max.counts <- counts[, sapply(.SD, max)]
pdf("pdf/cutoff_motif_density.pdf", 4, 4)
plot(density(counts[, sapply(.SD, max)]))
abline(v= 20)
dev.off()
keep <- names(max.counts)[max.counts<20]
counts <- counts[, keep, with= FALSE]
enr <- enr[name %in% keep]
enr[, name:= droplevels(name)]
enr[, motif:= droplevels(motif)]

# Add predicted activity to seq.info ----
seq.info[pred, pred.act:= i.vista_sigmoid, on= "id==vista_ID"]

# Add annotations for the different heatmaps ----
# First heatmap -> subset sequences
seq.info[, hm1:= seq(.N)<=200, label]

# Cluster ledidi sequences and select best clusters ----
# Clustering
seq.info[label %in% c("vista_ts", best.ledidi), annot.1:= ifelse(selected,
                                                                 paste0(best.ledidi, "_selected"),
                                                                 as.character(label))]
hm2 <- counts[seq.info$label %in% c("vista_ts", best.ledidi)]
set.seed(2)
pseudo <- min(unlist(hm2)[unlist(hm2)>0])/2
norm1 <- apply(hm2, 2, function(x) log2(x+pseudo))
k1 <- kmeans(norm1, centers = 9)

# Define color breaks
br <- seq(0, 3, length.out= 21)
Cc <- c("grey", colorRampPalette(c("pink", "red"))(length(br)-1))
leg <- "Motif/kb"

# Plot
vl_par(mai= c(.9,.5,.9,1.5))
vl_heatmap(x = hm2,
           cluster.rows = k1$cluster,
           row.annotations = as.character(na.omit(seq.info$annot.1)),
           breaks = br,
           col= Cc,
           legend.title = leg,
           show.rownames = FALSE)
seq.info[label %in% c("vista_ts", best.ledidi), k1:= k1$cluster]

# Select top sequences ----
seq.info[, annot.3:= fcase(!is.na(k1) & selected, paste0(best.ledidi, "_selected"),
                           !is.na(k1) & pred.act>=predict.act.cutoff, paste0("vista_ts>=", predict.act.cutoff),
                           default= NA_character_)]
# For each cluster, compute % of sequences containing each motif
perc <- counts[, lapply(.SD, function(x) sum(x>0)/length(x)), seq.info[, .(k1)]]
perc <- perc[!is.na(k1)]
perc <- melt(perc,
             id.vars = "k1",
             value.name = "perc",
             variable.name = "motif")
perc <- perc[perc>0.6, !"perc"]
# Get motif counts
top <- cbind(seq.info[, .(id, k1, annot.3)], counts)
top <- top[!is.na(k1) & !is.na(annot.3)]
top <- melt(top,
            id.vars = c("id", "k1", "annot.3"),
            value.name = "count",
            variable.name = "motif")
top <- merge(top,
             perc,
             by= c("k1", "motif"))
# Save motif enriched in designed enhancers
enh.mot <- unique(top[grepl(best.ledidi, annot.3), motif])
top[, mot.perc:= sum(count>0)/.N, id]
# Select top sequences containing the most of these motifs
setorderv(top, "mot.perc", -1)
top <- rbind(top[, .(id= unique(id)[1:5]), k1],
)
top <- rbind(top, data.table())
seq.info[, annot.3:= ifelse(id %in% top$id, annot.3, NA_character_)]
seq.info[, k3:= ifelse(!is.na(annot.3), k1, NA)]

# Function to plot barplot and heatmap 
barplots <- function(sub= enr$name) {
  # Contributions
  unique(enr[name %in% sub, .(name, enh.contrib, acc.contrib)])[, {
    # Enhancer contribution
    vl_barplot(enh.contrib,
               xaxt= "n",
               ylab= "Enh.contrib")
    # Accessibility contribution
    vl_barplot(acc.contrib,
               xaxt= "n",
               ylab= "Acc.contrib")
    
  }]
  # Plot barplots top heatmap
  log2OR <- dcast(enr[name %in% sub],
                  cl~name,
                  value.var = "log2OR")
  log2OR <- as.matrix(log2OR, 1)
  padj <- dcast(enr[name %in% sub],
                cl~name,
                value.var = "padj")
  padj <- as.matrix(padj, 1)
  # log2OR[padj>0.05] <- NA
  stars <- padj
  stars[padj<0.05 & log2OR>log2(1.5)] <- "*"
  stars[stars!="*"] <- NA
  # Plot enrichments
  vl_heatmap(log2OR,
             cluster.rows = F,
             show.numbers = stars,
             legend.title = "log2OR",
             breaks = seq(-3, 3, length.out= 21))
}

# Save final selected sequences ----
final <- list(seq= seq.info[annot.3==paste0(best.ledidi, "_selected"), .(id, width, cluster= k3)],
              cluster.motifs= factor(enh.mot, unique(enh.mot)))
saveRDS(final,
        "db/selected_sequences/heart_final_selection.rds")

# Heatmap individual sequences motif counts ----
pdf("pdf/motif_heatmap_selected_sequences.pdf", width = 14.5, height = 9)
vl_par(mai= c(.05, 1.25, .05, 1.6),
       omi= c(.5, 0, 0.1, 0),
       mgp= c(2.25, .55, 0.2),
       xaxs= "i",
       lwd= .5,
       cex.lab= .45)
layout(matrix(c(1,2,3,5,4)),
       heights = c(.5,.5,.3,.4,6))
# Heatmap motif counts
barplots()
vl_heatmap(counts[seq.info$hm1, ],
           cluster.rows = seq.info[(hm1), label],
           show.rownames = F,
           breaks = br,
           col= Cc,
           legend.title = leg,
           show.row.clusters = "left")
plot.new()
# Cluster vista and designed sequences
barplots()
vl_heatmap(counts[!is.na(seq.info$k1), ],
           cluster.rows = as.numeric(na.omit(seq.info$k1)),
           row.annotations = as.character(na.omit(seq.info$annot.1)),
           show.rownames = F,
           breaks = br,
           col= Cc,
           legend.title = leg,
           show.row.clusters = "left")
plot.new()
# Cluster vista and designed sequences
par(mai= c(.05, 3, .05, 3))
barplots(sub= unique(perc$motif))
par(mai= c(2.75, 3, .05, 3))
sub <- counts[!is.na(seq.info$k3), unique(perc$motif), with= FALSE]
sub <- as.matrix(sub)
rownames(sub) <- seq.info[!is.na(k3), id]
pl <- vl_heatmap(sub,
                 cluster.rows = as.numeric(na.omit(seq.info$k3)),
                 row.annotations = as.character(na.omit(seq.info$annot.3)),
                 show.rownames = TRUE,
                 breaks = br,
                 col= Cc,
                 legend.title = leg,
                 show.colnames = FALSE)
tiltAxis(pl$cols$x.pos,
         labels = colnames(sub),
         col= ifelse(colnames(sub) %in% enh.mot, "red", "black"))
pl$rows[, cluster:= as.integer(cluster)]
perc[, {
  x <- range(pl$cols[.BY, x.pos, on= "name==motif"])
  y <- range(pl$rows[.BY, y.pos, on= "cluster==k1"])
  rect(xleft = x[1]-.5,
       ybottom = y[1]-.6,
       xright = x[2]+.5,
       ytop= y[2]+.6,
       lwd= 1,
       lty= ifelse(motif %in% enh.mot, 1, 3),
       lend= 2,
       xpd= NA,
       border= "black")
}, .(motif, k1)]
dev.off()

