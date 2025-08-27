setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import sequence info ----
dat <- list(heart= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_heart.rds"),
            limb= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_limb.rds"),
            midbrain= readRDS("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_midbrain.rds"))
dat <- rbindlist(dat, idcol = "tissue")


# Set cutoffs ----
cutoffs <- c(heart= 6, midbrain= 7, limb= 5)

# ledidi parameters ----
best.ledidi <- "ledidi_12_14"

# Plotting parameters ----
Cc <- c("red", "limegreen", "darkgrey", "cornflowerblue")

# Density predictions ----
pdf("pdf/density_predictions_all_tissues_evgeny.pdf", width = 7/3*5, height = 6.5)
vl_par(mfrow= c(3, 5),
       omi= c(0,0,.5,0))
for(tiss in c("heart", "limb", "midbrain")) {
  # Compute tissues stats
  ledidi <- dat[tissue==tiss & label==best.ledidi]
  
  # Plot blast pie chart
  bl <- table(ifelse(ledidi$blast, "BLAST hit", "No BLAST hit"))
  par(mai= c(.6, .6, .6, .6))
  pie(bl,
      labels = paste0(names(bl), "\n(", bl, ")"),
      main= paste("BLAST", tiss),
      cex= .8)
  
  # For each model ----
  for(model in c("heart", "limb", "midbrain")) {
    
    # Ledidi designed activities
    d1 <- density(ledidi[[paste0("act_", model)]], na.rm= T)
    
    # Vista tiles within tissue
    if(model == tiss) {
      # Ledidi designed
      d2 <- density(dat[tissue==tiss & label=="vista_ts"][[paste0("act_", model)]], na.rm= T)
      # Inactive
      d3 <- density(dat[tissue==tiss & label=="vista_inactive_all_tissues"][[paste0("act_", model)]], na.rm= T)
      # X and Y lim
      xlim <- range(c(d1$x, d2$x, d3$x))
      ylim <- range(c(d1$y, d2$y, d3$y))
      
    } else {
      
      # Vista ts other tissues
      d4 <- density(dat[tissue==model & label=="vista_ts"][[paste0("act_", model)]], na.rm= T)
      xlim <- range(c(d1$x, d4$x))
      ylim <- range(c(d1$y, d4$y))
      
    }
    
    # Initiate plot ----
    par(mai= c(.5, .3, .2, .6))
    plot(d1,
         xlim= xlim,
         ylim= c(ylim[1], ylim[2]+diff(ylim)/5), # Add space labels
         main= if(tiss=="heart") paste("Model:", model) else NA, # 1st line
         col= "red",
         lwd= 2,
         ylab= "Density",
         xlab= NA,
         xaxt= "n")
    
    # Axes labels ----
    axis(1, padj= -1.25)
    title(xlab= "Predicted activity", line = .75)
    if(model== "limb") # Middle plot
      mtext(text = paste0(tiss, " sequences"), cex= 1.2, line = 1.5)
    
    # Add other densities ----
    if(tiss==model) {
      # Density
      lines(d2, col= "limegreen", lwd= 2)
      lines(d3, col= "darkgrey", lwd= 2)
      abline(v= cutoffs[tiss], lty= 3)
      # Legend
      leg <- c(paste0(tiss, c(" ledidi", " vista ts")), "vista inact.")
      vl_legend(x.adj = -2,
                col= c("red", "limegreen", "darkgrey"),
                lwd= 2,
                legend = leg)
      # Activity cutoff
      text(cutoffs[tiss],
           c(par("usr")[4]-c(1,2)*strheight("M")),
           c(paste0(">", cutoffs[tiss]), paste0(sum(ledidi$active), "/1200")),
           pos= 4,
           cex= .8)
    } else {
      
      # Density
      lines(d4, col= "limegreen", lwd= 2)
      leg <- c(paste0(tiss, " ledidi"), paste0(model, " vista ts"))
      abline(v= 0, lty= 3)
      # Legend
      vl_legend(x.adj = -2,
                col= c("red", "limegreen"),
                lwd= 2,
                legend = leg)
      # Specificity cutoff
      text(0,
           c(par("usr")[4]-c(1,2)*strheight("M")),
           c("<0", paste0(sum(ledidi[, paste0("act_", model), with= FALSE]<0), "/1200")),
           pos= 2,
           cex= .8)
    }
  }
  
  # Compute final pie chart ----
  summary <- ledidi[, fcase(blast & !active & !specific, "blast+inact+aspecific",
                            blast & !active, "blast+inact",
                            blast & !specific, "blast+aspecific",
                            !active & !specific, "inact+aspecific",
                            !active, "inactive",
                            !specific, "aspecific",
                            blast, "blast",
                            default= "Passed filtering")]
  summary <- factor(
    summary,
    c("Passed filtering", "blast", "aspecific",
      "inactive", "blast+aspecific", "blast+inact",
      "inact+aspecific", "blast+inact+aspecific")
  )
  summary <- table(summary)
  par(mai= c(.6, .4, .6, .8))
  Cc <- adjustcolor(rainbow(length(summary)), .6)
  
  # Plot ----
  pie(summary,
      labels = ifelse(summary>50, paste0(names(summary), "\n(", summary, ")"), NA),
      main= paste("Summary", tiss),
      cex= .6,
      col= Cc)
  vl_legend(fill= Cc,
            legend = names(summary),
            cex= .5,
            x.adj = .5)
  print(summary)
}
dev.off()