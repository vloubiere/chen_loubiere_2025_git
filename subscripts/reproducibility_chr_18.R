setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)
source("git_deepATAC/function/compute_AUC.R")
source("git_deepATAC/function/compute_PPV.R")

# Import metadata ----
meta <- readRDS("Rdata/metadata_ATACSeq_models_processed.rds")

# Retrieve prediction files ----
ATACfolder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/model/ATAC_model/"
meta[, pred_ATAC_test:= file.path(ATACfolder, ID, tissue, paste0("results_", fold, "_", tissue, "_DeepSTARR2_", replicate), paste0(fold, "_sequences_test.fa_predictions_enhancer_Model.txt")), .(ID, tissue, fold, replicate)]
meta[, pred_ATAC_chr18:= file.path(ATACfolder, ID, tissue, paste0("results_", fold, "_", tissue, "_DeepSTARR2_", replicate), "whole_chr18.fa_predictions_enhancer_Model.txt"), .(ID, tissue, fold, replicate)]
VISTAfolder <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/model/VISTA_model/"
meta[, pred_VISTA_test:= file.path(VISTAfolder, ID, tissue, paste0("results_", fold, "_", tissue, "_DeepSTARR2_", replicate), paste0(fold, "_sequences_test.fa_predictions_enhancer_Model.txt")), .(ID, tissue, fold, replicate)]
meta <- meta[file.exists(pred_ATAC_test) & file.exists(pred_ATAC_chr18) & file.exists(pred_VISTA_test)]

# Import predicted activities ----
observed <- meta[ID==ID[1], {
  .c <- lapply(obs_ATAC_test, fread)
  names(.c) <- paste0(tissue, "_", fold, "_", replicate)
  .c <- rbindlist(.c, idcol = T)
  .c <- .c[grepl("chr18", ID)]
  .c <- dcast(.c, ID~.id, value.var = "score")
  .c <- na.omit(.c[, -1])
}, .(ID)]

# Compare PCC observed vs predicted for each model
Cc <- c("black", "midnightblue", viridis::viridis(9), "tomato", "red")
printFUN <- function(x) apply(x, c(1,2), function(x) gsub("^0", "", as.character(round(x, 2))))

# Plot pcc
pdf("pdf/ATAC_gradient_reproducibility_chr18.pdf", 30, 26)
vl_par(mai= c(1.1,.9,.25,.75), mfrow= c(5,5), lwd= .25)
hm <- vl_heatmap(cor(observed[,-1]),
                 main= "Observed",
                 tilt.colnames = T,
                 legend.cex = .8,
                 display.numbers = T,
                 display.numbers.cex = .8,
                 display.numbers.FUN = printFUN,
                 col= Cc,
                 breaks = seq(.7, 1, length.out= length(Cc)))
for(i in 1:4)
  plot.new()
vl_par(mai= c(.6,.4,.25,.75), cex= .66)
meta[, {
  .c <- lapply(pred_ATAC_test, fread)
  names(.c) <- paste0(tissue, "_", fold, "_", replicate)
  .c <- rbindlist(.c, idcol = T)
  .c <- .c[grepl("chr18", location)]
  .c <- dcast(.c, location~.id, value.var = "Predictions")
  .c <- na.omit(.c)
  hm$x <- cor(.c[,-1])
  plot(hm,
       show.rownames = FALSE,
       show.colnames = FALSE,
       show.row.dendrogram = FALSE,
       show.col.dendrogram = FALSE,
       legend.cex = .8,
       display.numbers = T,
       display.numbers.cex = .8,
       display.numbers.FUN = function(x) apply(x, c(1,2), function(x) gsub("^0", "", as.character(round(x, 2)))),
       main= gsub("^model1_bulkATAC_", "", ID),
       col= Cc,
       breaks = seq(.7, 1, length.out= length(Cc)))
  print(".")
}, ID]
dev.off()

meta[, {
  .c <- lapply(pred_ATAC_test, fread)
  names(.c) <- paste0(tissue, "_", fold, "_", replicate)
  .c <- rbindlist(.c, idcol = T)
  .c <- .c[grepl("chr18", location)]
  .c <- dcast(.c, location~.id, value.var = "Predictions")
  .c <- na.omit(.c)
  vl_heatmap(cor(.c[,-1]))
  print(".")
}, .(ID, augmentation, balancing, weight)]

# Observed activities ----
meta[, c("obs", "class"):= {
  ATAC <- fread(obs_ATAC_test)
  VISTA <- fread(obs_VISTA_test)
  .(.(c(ATAC$score, # ATAC fold-specific test
        log2(obsChr18[[tissue]][ok]+1), # ATAC Chr18
        VISTA$score)), # VISTA
    .(c(ATAC$class,
        rep("chr18", length(ok)),
        VISTA$class)))
}, .(tissue, obs_ATAC_test, obs_VISTA_test)]

# Predicted activities
meta[, pred:= .(
  .(c(fread(pred_ATAC_test)$Predictions, # ATAC fold-specific test
      fread(pred_ATAC_chr18)$Predictions[ok], # ATAC Chr18
      fread(pred_VISTA_test)$Predictions)) # VISTA
), .(pred_ATAC_test, pred_ATAC_chr18, pred_VISTA_test)]

# Compute reproducibility metrics ----
stats <- meta[, {
  # Vars
  obs <- obs[[1]]
  pred <- pred[[1]]
  cl <- class[[1]]
  
  # Classes
  foldSpe <- cl %in% c("globallyClosed", "globallyOpen", "specificClosed", "tissueSpecific") # fold-specific ATAC test set
  tissueSpe <- cl %in% c("specificClosed", "tissueSpecific") # Tissue-spe ATAC test
  chr18 <- cl == "chr18" # Whole chromosome 18 ATAC
  vista <- cl %in% c("active", "inactive")
  
  # Compute PPV
  PPV <- plot_predictive_positive(label = obs[vista],
                                  predicted = pred[vista])
  
  # Return metrics
  .(pccGlobal= cor(obs[foldSpe], pred[foldSpe]), # PCC fold-specific test set
    pccSpecific= cor(obs[tissueSpe], pred[tissueSpe]), # PCC tissue-specific peaks
    pccChr18= cor(obs[chr18], pred[chr18]), # PCC whole chromosome 18
    rocAUC= vl_ROC_AUC(label= obs[vista], predicted = pred[vista]), # rocAUC VISTA tiles
    prAUC= PRROC::pr.curve(scores.class0 = pred[vista][obs[vista]==1], # prAUC VISTA tiles
                           scores.class1 = pred[vista][obs[vista]==0],
                           curve = TRUE)$auc.integral,
    min_PPV= PPV$min_PPV,
    max_PPV= PPV$PPV_at_cutoff,
    delta_PPV= PPV$PPV_at_cutoff/PPV$min_PPV)
}, .(tissue, fold, replicate, augmentation, balancing, weight)]

# Melt ----
stats <- melt(stats,
              id.vars = c("tissue", "fold", "replicate", "augmentation", "balancing", "weight"),
              measure.vars = c("pccGlobal", "pccSpecific", "pccChr18", "rocAUC", "prAUC", "min_PPV", "max_PPV", "delta_PPV"))
stats[, weight:= factor(weight, c("no", "equal", "TS", "forced"))]
setorderv(stats,
          c("augmentation", "balancing", "weight"))
stats[, weight:= as.character(weight)]

# Plotting function ----
pl <- substitute(
  {
    .c <- .SD[, .(mean(value), .(value), name= paste(unlist(.BY), collapse= " ")), .(augmentation, balancing, weight)]
    xaxt <- ifelse(variable=="delta_PPV", "s", "n")
    yaxt <- ifelse(tissue=="forebrain", "s", "n")
    ylim <- if(variable=="delta_PPV")
      c(0, 10) else if(grepl("PPV", variable))
        c(0, 100) else c(0, 1)
    roundFUN <- if(grepl("PPV", variable))
      function(x) round(x, 1) else
        function(x) gsub("^0", "", as.character(round(x, 3)))
    .c[, {
      vl_barplot(V1,
                 bar.labels = roundFUN(V1),
                 bar.labels.cex = .5,
                 individual.var = V2,
                 ylim= ylim,
                 xaxt= xaxt,
                 yaxt= yaxt,
                 names.arg = name,
                 ind.col = rep(c("tomato", "limegreen", "cornflowerblue"), each= 2),
                 space= c(0.1, rep(c(.1,.1,.1,.5), 4))[seq(V1)])
      if(tissue=="forebrain")
        title(ylab= variable, xpd= NA)
      if(variable=="pccGlobal")
        title(main= tissue[1], xpd= NA, line = 2)
    }]
  }
)

# PDF ----
pdf("pdf/ATAC_gradient_reproducibility_metrics.pdf", width = 25, height= 18)
mat <- matrix(1:(6*8), ncol= 6, byrow = T)
layout(mat)
vl_par(mai= c(.1,.2,.2,0),
       omi= c(1,1,1,0))
stats[, {
  # Global PCC ----
  eval(pl)
}, .(tissue, variable)]
dev.off()

# Average replicates ----
dat <- meta[, {
  .(obs= .(obs[[1]]),
    class= .(class[[1]]),
    pred= .(rowMeans(do.call(cbind, pred))))
}, .(ID, tissue, fold)]

# Group all folds ----
dat <- dat[, {
  .(obs= .(unlist(obs)),
    class= .(unlist(class)),
    pred= .(unlist(pred)))
}, .(ID, tissue)]

# Save data ----
saveRDS(dat,
        "db/predictions/model1_20240902_ATACSeq_gradient.rds")
