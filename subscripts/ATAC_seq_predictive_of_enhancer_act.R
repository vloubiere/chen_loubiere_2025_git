setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metadata ----
meta <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/metadata_scATACSeq_model_predictions.rds")

# Import labels ----
labels <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/augmentation/snATAC/vistaTiles_1001_50.rds")

# Add accessibility ----
acc <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/scATAC_with_signals.rds")
labels[acc, acc:= i.limb, on= c("seqnames", "start<=end", "end>=start")]
labels[is.na(acc), acc:= 0]

# Import data ----
dat <- meta[tissue=="limb" & dataset=="snATACseq" & ID=="model1_snATACseq_tsx3_nogcAug_1xBal_noW", fread(pred_scATAC_vista), pred_scATAC_vista]
dat <- dat[, .(pred= mean(Predictions)), location]
# Remove human sequences
dat <- dat[!grepl("^hs", location)]
dat[labels, enh:= i.limb, on= "location==ID"]
dat[labels, obs:= i.acc, on= "location==ID"]

dat[, class:= fcase(enh==1 & obs==1, "Open enhancer",
                    enh==1 & obs==0, "Closed enhancer",
                    enh==0 & obs==1, "Open inactive",
                    default = "Closed inactive")]
dat[, class:= factor(class, c("Open enhancer",
                              "Open inactive",
                              "Closed enhancer",
                              "Closed inactive"))]

Cc <- c("red", "blue", "green", "grey")
plot(c(0,1.8), c(0,10), type= "n")
dat[, {
  lines(density(pred), col= Cc[class])
}, class]
legend("topright",
       fill= Cc,
       legend = levels(dat$class))
