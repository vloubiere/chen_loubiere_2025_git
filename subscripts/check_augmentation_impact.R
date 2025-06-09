setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
source("git_deepATAC/function/compute_AUC.R")

# Import data ----
x <- readRDS("/groups/stark/shenzhi.chen/model/deepstarr2-mouse/Data/Input/single_task/Muti_tissue_accessibility_model_20240617/db/VISTA_20240617_test_performance_enhancer.rds")

# Variables ----
predScores <- c("Forebrain",
                "Heart",
                "Hindbrain",
                "Limb",
                "Midbrain",
                "Neural_tube")
tissues <- tolower(predScores)
predScores <- paste0(rep(predScores, each= 3),
                     c("_aug_1", "_aug_2", "_aug_3"))

# Extract labels
AUC <- merge(x[, .(ID2, fold)],
             melt(x,
                  id.vars = c("ID2", "fold"),
                  measure.vars = tissues,
                  variable.name = "tissue",
                  value.name = "label"),
             by= c("ID2", "fold"),
             all.y= T)
# Extract scores
scores <- melt(x, 
               id.vars = c("ID2", "fold"),
               measure.vars = predScores,
               variable.name = "replicate",
               value.name = "score")
scores[, c("tissue", "replicate"):= .(tolower(gsub("_[^_]*$", "", replicate)), gsub(".*(.{1})$", "\\1", replicate))]
AUC <- merge(AUC,
             scores,
             by= c("ID2", "fold", "tissue"),
             all.y= T)
AUC[, label:= as.logical(label=="Active")]
# Order
setorderv(AUC,
          c("tissue", "fold", "replicate", "score"),
          c(1, 1, 1, -1))