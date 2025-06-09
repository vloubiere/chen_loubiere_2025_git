setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Save bed file (originally in rds) ----
bed <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/augmentation/snATAC/vistaTiles_1001_50.rds")
tmp <- tempfile(fileext = ".bed")
exportBed(bed[, .(seqnames, start, end, name= ID)], tmp)

# Import contributions scores ----
contrib <- importContrib(h5 = "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/model/VISTA_model/model1_bulkATAC_tsx3Aug_2xBal_noW/heart/results_fold01_heart_DeepSTARR2_rep1/Model_all_vista_seq.fa_dinuc_shuffle_deepSHAP_DeepExplainer_importance_scores.h5",
                         bed = tmp,
                         fa = "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/training_dataset/VISTA_snATAC_20241001/all_vista_seq.fa",
                         FUN = function(x) colSums(x))
# npz_file <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/model/VISTA_model/model1_bulkATAC_tsx3Aug_2xBal_noW/heart/results_fold01_heart_DeepSTARR2_rep1/all_vista_seq_contrib.npz"

# Simplify sequences ----
uniq <- contrib[grepl("_0$|_1000$|_2000$|_3000$|_4000$|_5000$", name)]

# Call seqlets ----
seqLets <- vlite::contribSeqlets(uniq)
enr <- seqLets[log2OR>1 & FDR<1e-5]
sel <- enr[, .N, seqlvls][order(-N)]$seqlvls[2]

vlite::contribSeqLogo(uniq[name==sel],
                      start = 600,
                      end= 800,
                      highlight = seqLets[seqlvls==sel, seq(start, end), .(start, end)]$V1)
