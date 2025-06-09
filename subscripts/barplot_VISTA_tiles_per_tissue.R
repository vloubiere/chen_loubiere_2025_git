setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import usable ones ----
dat <- readRDS("Rdata/folds_training_valid_test_sets_cov_act.rds")[modelType=="VISTA"]

# Count enhancers ----
count <- melt(dat,
              id.vars = "ID",
              measure.vars = c("heart", "limb", "forebrain", "hindbrain", "midbrain", "neuralTube"))
count[, vistaID:= tstrsplit(ID, "[.]", keep= 1)]
count <- unique(count[, !"ID"])
count[, class:= ifelse(sum(value)==0, "GloballyInactive", as.character(NA)), vistaID]
count[, class:= ifelse(all(value==1), "GloballyActive", class), vistaID]
count[value==1 & is.na(class), class:= "Active"]
count[value==0 & is.na(class), class:= "Inactive"]
count[, class:= factor(class, 
                       c("Active", "Inactive", "GloballyActive", "GloballyInactive"))]

mat <- dcast(count, class~variable)
mat <- as.matrix(mat, 1)


Cc <- c("limegreen", "cornflowerblue", "red", "grey")
pdf("pdf/VISTA_N_enhancers_per_tissue.pdf", 3.3, 3)
vl_par(las= 2,
       mai= c(.9, .7, .9, 1.2))
barplot(mat,
        col= adjustcolor(Cc, .3),
        ylab= "N enhancers",
        border= NA)
legend(par("usr")[2],
       par("usr")[4],
       fill= adjustcolor(Cc, .3),
       legend = rownames(mat),
       xpd= NA,
       bty= "n",
       cex= .7,
       border= NA)
dev.off()