setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import sequence info ----
dat <- list(heart= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_heart.rds"),
            limb= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_limb.rds"),
            midbrain= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_midbrain.rds"))
dat <- rbindlist(dat, idcol = "tissue")


# Plot activity densities ----
Cc <- c("red", "limegreen", "darkgrey", "cornflowerblue")


pdf("pdf/density_predictions_all_tissues.pdf", width = 3.5, 6)
vl_par(mfrow= c(3, 1),
       mai= c(.5, .5, .2, 1.7),
       omi= c(0,0,.5,0))
for(tiss in c("heart", "limb", "midbrain")) {
  for(model in c("heart", "limb", "midbrain")) {
    # Vista TS tissue
    d1 <- density(dat[tissue==tiss & label=="ledidi_12_14"][[paste0("act_", model)]], na.rm= T)
    # Vista ts other tissue
    if(model == tiss) {
      # Ledidi designed
      d2 <- density(dat[tissue==tiss & label=="vista_ts"][[paste0("act_", model)]], na.rm= T)
      # Inactive
      d3 <- density(dat[tissue==tiss & label=="vista_inactive_all_tissues"][[paste0("act_", model)]], na.rm= T)
      # X and Y lim
      xlim <- range(c(d1$x, d2$x, d3$x))
      ylim <- range(c(d1$y, d2$y, d3$y))
    } else {
      # vista other tissues
      d4 <- density(dat[tissue==model & label=="vista_ts"][[paste0("act_", model)]], na.rm= T)
      xlim <- range(c(d1$x, d4$x))
      ylim <- range(c(d1$y, d4$y))
    }
    ylim <- c(ylim[1], ylim[2]+diff(ylim)/10)
    # Plot
    plot(d1,
         xlim= xlim,
         ylim= ylim,
         main= paste("Model:", model),
         col= "red",
         lwd= 2,
         ylab= "Density",
         xlab= "Predicted activity")
    mtext(text = paste0(tiss, " sequences"), outer = T, cex= 1.5, line = 1.5)
    if(tiss==model) {
      # Density
      lines(d2, col= "limegreen", lwd= 2)
      lines(d3, col= "darkgrey", lwd= 2)
      cutoff <- switch(tiss,
                       "heart"= 6,
                       "midbrain"= 7,
                       "limb"= 5)
      abline(v= cutoff, lty= 3)
      # Legend
      leg <- paste0(tiss, c(" ledidi", " vista ts", "vista inact."))
      vl_legend(col= c("red", "limegreen", "darkgrey"),
                lwd= 2,
                legend = leg)
      text(cutoff,
           c(par("usr")[4]-c(1,2)*strheight("M")),
           c(paste0(">", cutoff), paste0(sum(dat[label=="ledidi_12_14" & tissue==tiss, selected]), "/1200")),
           pos= 4,
           cex= .8)
    } else {
      # Density
      lines(d4, col= "limegreen", lwd= 2)
      leg <- c(paste0(tiss, " ledidi"), paste0(model, " vista ts"))
      abline(v= 0, lty= 3)
      # Legend
      vl_legend(col= c("red", "limegreen"),
                lwd= 2,
                legend = leg)
      text(0,
           par("usr")[4]-strheight("M"),
           "<0",
           pos= 2,
           cex= .8)
    }
  }
}
dev.off()

# dat[, selected:= (tissue=="heart" & act_heart>7 & act_limb<(0) & act_midbrain<(0) & !blast)
#     | (tissue=="limb" & act_heart<(0) & act_limb>5 & act_midbrain<(0) & !blast)
#     | (tissue=="midbrain" & act_heart<(-1) & act_limb<(-1) & act_midbrain>7 & !blast)]
# dat[, pred.act:= get(paste0("act_", tissue)), tissue]
# dat[, sum(selected, na.rm = T), tissue]
# 
# vl_par()
# vl_boxplot(pred.act~selected+tissue,
#            dat,
#            tilt.names= T,
#            ylab= "Predicted activity")
