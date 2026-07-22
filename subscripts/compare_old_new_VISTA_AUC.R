setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
source("git_deepATAC/function/compute_AUC.R")

# Import data ----
old <- readRDS("/groups/stark/shenzhi.chen/model/deepstarr2-mouse/Data/Input/single_task/Muti_tissue_accessibility_model_20240617/db/VISTA_old_test_performance_enhancer.rds")
new <- readRDS("/groups/stark/shenzhi.chen/model/deepstarr2-mouse/Data/Input/single_task/Muti_tissue_accessibility_model_20240617/db/VISTA_20240617_test_performance_enhancer.rds")

# Compute AUC for each fold/model/tissue ----
oldAUC <- computeAUC(old)
newAUC <- computeAUC(new)
comp <- list(old= oldAUC,
             new= newAUC)
comp <- rbindlist(comp, idcol = T)

# Plot ----
pdf("pdf/AUC_old_vs_new_VISTA_db.pdf", 3, width = 16)
vl_par(mfrow= c(2.5, 11),
       mai= c(.3,.3,.2,.2),
       omi= c(0,1,.5,0))
comp[, {
  .SD[, {
    .SD[, {
      plot(c(0, 1),
           c(0,1),
           type= "n",
           xlab= "FPR",
           ylab= "TPR")
      if(.GRP==1)
        text(par("usr")[1]-strwidth("MMM"),
             mean(par("usr")[c(3,4)]),
             pos= 2,
             .id,
             xpd= NA,
             cex= 2)
      .SD[, {
        lines(FPR[[1]],
              TPR[[1]],
              col= "red")
        abline(0,
               1,
               lty= 3)
      }, replicate]
    }, fold]
    vl_boxplot(AUC~fold,
               tilt.names= T,
               ylab= "AUC",
               ylim= c(0.4, .8))
    abline(h= .5, lty= 3)
    .SD
  }, .id]
  mtext(tissue, outer = T)
  print(".")
}, tissue]
dev.off()