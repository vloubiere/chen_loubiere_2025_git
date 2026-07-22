setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import metadata ----
dat <- readRDS("db/predictions/20240807_ATACSeq_model_testSet_predictions.rds")
dat <- dat[fold %in% paste0("fold0", 1:7)]

# Cast
obsHm <- dcast(dat, ID~tissue, value.var = "observed", fun.aggregate = function(x) mean(x, na.rm =T))
predHm <- dcast(dat, ID~tissue, value.var = "predicted", fun.aggregate = function(x) mean(x, na.rm =T))
Hm <- merge(obsHm,
            predHm,
            by= "ID",
            suffixes= c(".obs", ".pred"))
labels <- readRDS("Rdata/folds_training_valid_test_sets_cov_act.rds")
Hm[labels, class:= i.class, on= "ID"]
Hm <- na.omit(Hm)

# # Scale data ---
# # Assuming table1 is your reference table
# means <- apply(Hm[, forebrain.obs:neuralTube.obs], 2, mean, na.rm = TRUE)
# sds <- apply(Hm[, forebrain.obs:neuralTube.obs], 2, sd, na.rm = TRUE)
# 
# # Assuming table2 is the table you want to scale
# scaledObs <- scale(Hm[, forebrain.obs:neuralTube.obs], center = means, scale = sds)
# scaledPred <- scale(Hm[, forebrain.pred:neuralTube.pred], center = means, scale = sds)

# Plot ----
pdf("pdf/ATAC_heatmap_obs_exp.pdf",
    height = 5.5,
    width = 11.5)
vl_par(mai= c(1.2, 2.5,.9,2),
       mfrow= c(1, 2),
       cex.axis= 1)
vl_heatmap(scale(Hm[, forebrain.obs:neuralTube.obs]),
           row.clusters= Hm$class,
           cluster.cols = F,
           show.rownames = F,
           breaks = c(-1.5, 0, 1.5),
           # col= c("blue", "yellow"),
           kmeans.k= 7,
           row.clusters.pos = "left",
           tilt.colnames = T,
           show.legend = F,
           main= "Observed")
vl_heatmap(scale(Hm[, forebrain.pred:neuralTube.pred]),
           row.clusters= Hm$class,
           cluster.cols = F,
           show.rownames = F,
           breaks = c(-1.5, 0, 1.5),
           # col= c("blue", "yellow"),
           kmeans.k= 7,
           row.clusters.pos = "left",
           tilt.colnames = T,
           legend.title = "ATAC-Seq Z-score",
           main= "Predicted")
dev.off()