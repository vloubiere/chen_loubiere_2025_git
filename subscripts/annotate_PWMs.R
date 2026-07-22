setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import motifs ----
mot <- list(heart= readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/final_enhancer_selection/all_sequences_heart_contri_score.rds"),
              limb= readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/final_enhancer_selection/all_sequences_limb_contri_score.rds"),
              midbrain= readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/final_enhancer_selection/all_sequences_midbrain_contri_score.rds"))
mot <- rbindlist(mot, fill = T, idcol = "tissue")
mot[, c("limb", "midbrain", "heart"):= NULL]
setnames(mot, c("contrib", "access"), c("enh.contrib", "acc.contrib"))
mot <- dcast(mot, motif~tissue, value.var = list("enh.contrib", "acc.contrib"))

# Add cluster ----
clust <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/final_enhancer_selection/all_motifs_clusters_metatable.rds")
mot[clust, cluster:= i.cluster, on= "motif"]
setcolorder(mot, "cluster")

# Rename relevant motif clusters ----
# Accessibility
mot[grepl("CREB", cluster, ignore.case = T), new.name:= paste0("CREB.accessibility")]
mot[grepl("AP1", cluster, ignore.case = T), new.name:= paste0("AP1.accessibility")]
mot[grepl("ZIC", cluster, ignore.case = T), new.name:= paste0("ZIC.accessibility")]
mot[grepl("KLF/SP/", cluster, ignore.case = T), new.name:= paste0("KLF/SP.accessibility")]
# Heart
mot[grepl("NKX2.5", motif, ignore.case = T), new.name:= "NKX2-5.heart"]
mot[grepl("HOXB2|HOXB3|HOXB5", motif, ignore.case = T), new.name:= "HOXB2-5.heart"]
mot[grepl("IRX", motif, ignore.case = T), new.name:= "IRX.heart"]
mot[grepl("TBX4", motif, ignore.case = T), new.name:= "TBX4.heart"]
mot[cluster %in% c("TFAP2/1", "TFAP2/2"), new.name:= "TFAP2.heart"]
mot[cluster %in% c("GATA", "MEF2", "TEAD", "EGR", "NRF1", "HAND1", "Ebox/CACGTG/1", "Ebox/CACGTG/2"),
      new.name:= paste0(cluster, ".heart")]
mot[cluster %in% c("PAX-halfsite", "PRDM16"), new.name:= paste0(cluster, ".unexpected")]
# Limb
mot[grepl("HOXD9|HOXD10|HOXD11|HOXD13", motif, ignore.case = T), new.name:= "HOXD9-13.limb"]
mot[grepl("GLI3", motif, ignore.case = T), new.name:= "GLI3.limb"]
mot[grepl("LMX1B", motif, ignore.case = T), new.name:= "LMX1B.limb"]
mot[grepl("SOX9", motif, ignore.case = T), new.name:= "SOX9.limb"]
mot[cluster %in% c("RUNX/1","RUNX/2"), new.name:= "RUNX.limb"]
mot[cluster %in% c("Ebox/CAGCTG", "Ebox/CAGATGG"), new.name:= "TWIST1.limb"]
mot[grepl("olig3", motif, ignore.case = T), new.name:= "olig3.limb"]
mot[grepl("HTF4", motif, ignore.case = T), new.name:= "HTF4.limb"]
# Midbrain
mot[grepl("OTX2", motif, ignore.case = T), new.name:= "OTX2.midbrain"]
mot[grepl("SOX9", motif, ignore.case = T), new.name:= "SOX9.midbrain"]
mot[grepl("OTX1", motif, ignore.case = T), new.name:= "OTX1.midbrain"]
mot[grepl("RFX", cluster, ignore.case = T), new.name:= "RFX.midbrain"]
mot[grepl("POU/1", cluster, ignore.case = T), new.name:= "POU.midbrain"]
mot[grepl("nkx21", motif, ignore.case = T), new.name:= "NKX21.midbrain"]
mot[grepl("NEUROD1", motif, ignore.case = T), new.name:= "NEUROD1.midbrain"]
mot[grepl("OLIG2", motif, ignore.case = T), new.name:= "OLIG2.midbrain"]
mot[grepl("Pax6", motif, ignore.case = T), new.name:= "PAX6.midbrain"]
mot[grepl("FOXG1", motif, ignore.case = T), new.name:= "FOXG1.midbrain"]
mot[grepl("Ngn2", motif, ignore.case = T), new.name:= "NGN2.midbrain"]
mot[grepl("TBR1_TBX_2", motif, ignore.case = T), new.name:= "TBR1.midbrain"]
mot[grepl("FEZF1", motif, ignore.case = T), new.name:= "FEZF.midbrain"]
mot[is.na(new.name), new.name:= paste0(cluster, ".other")]

# Clean ----
mot <- mot[!is.na(new.name)]
mot[, c("cluster", "annot"):= tstrsplit(new.name, "[.]", type.convert = T)]
mot$new.name <- NULL

# Add PWMs ----
pwm <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/jeff_non-redundant_pwms_jaspar_all.rds")
mot[, pwm:= as.list(pwm[match(mot$motif, sapply(pwm, TFBSTools::name))])]
pwm.perc <- readRDS("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/jeff_non-redundant_pwms_jaspar_all_percentage.rds")
mot[, pwm.perc:= as.list(pwm.perc[match(mot$motif, sapply(pwm.perc, TFBSTools::name))])]

# Save ----
saveRDS(mot,
        "Rdata/annotated_PWMs.rds")
