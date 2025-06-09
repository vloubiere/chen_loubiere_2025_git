setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Import data ----
dat <- readRDS("db/predictions/20240807_VISTA_model_testSet_predictions.rds")[aug=="No" & tissue=="heart"]
dat[, ID:= tstrsplit(ID, "_", keep= 1)]

# Retrieve sequence and compute CG content ----
seq <- readRDS("Rdata/folds_training_valid_test_sets_cov_act.rds")[modelType=="VISTA"]
seq[, nCG:= sapply(sequence, function(x) length(unlist(gregexpr("CG", x))))]

# Merge ----
dat[seq, nCG:= i.nCG, on= "ID"]

# Compute FDR ----
dat[, predPos:= predicted>0.7]
dat[, CGcut:= cut(nCG, quantile(nCG, seq(0, 1, length.out= 6)), include.lowest= T)]
res <- dat[, .(predPos= sum(predPos),
               V1= sum(predPos & (active))/sum(predPos)), .(fold, rep, CGcut)]
res <- res[predPos>20]

pdf("pdf/CG_content_PPV.pdf", 3, 3)
vl_par()
vl_boxplot(V1~CGcut,
           res,
           xlab= "Number of CG dinucleotides",
           ylab= "Positive Predicted Value",
           outline= T,
           tilt.names= T,
           main= "Heart (pos= pred>0.7)")
dev.off()
