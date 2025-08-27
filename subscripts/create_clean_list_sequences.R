setwd("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/")
# load library
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import 
for(tiss in c("heart","limb","midbrain")){
  
  # Import evo design and genomic sets
  act.acc <- readRDS(paste0("Rdata/motifs_enrichment_analysis/two_step_designed_sequences/all_sequences_for_motifs_analysis", tiss, ".rds"))
  act.acc <- split(act.acc[, .(sequence, id= names)], act.acc$label)
  
  # Import ledidi input sequences
  ledidi.input <- seqinr::read.fasta(
    paste0("result/sequence_design/Leididi_design/", tiss, "/all_sequences/all_ledidi_designed_sequence_ini.fasta"), as.string = T)
  ids <- gsub(".0_", "_", names(ledidi.input))
  ids <- unlist(tstrsplit(ids, paste0("_", tiss), keep= 1))
  ids <- gsub("designed_", "", ids)
  ledidi.input <- data.table(sequence= toupper(as.character(ledidi.input)),
                             id= rowid(ids))
  ledidi.input <- split(ledidi.input, paste0(ids, "_ini"))
  
  # Import ledidi designed sequences
  ledidi.design <- seqinr::read.fasta(
    paste0("result/sequence_design/Leididi_design/", tiss, "/all_sequences/all_ledidi_designed_sequence.fasta"),
    as.string = T)
  seq.names <- names(ledidi.design)
  ledidi.design <- data.table(sequence= toupper(as.character(ledidi.design)),
                              id= rowid(ids))
  # Add selected sequences
  ledidi.sel <- readRDS(
    paste0("Rdata/final_enhancer_selection/merged_all_selected_ledidi_designed_sequences_", tiss, ".rds"))
  ledidi.design[, selected:= seq.names %in% ledidi.sel]
  ledidi.design <- split(ledidi.design, ids)
  
  # Make final object ----
  obj <- c(
    list(rdm= act.acc[["rdm"]],
         vista_inactive_all_tissues= act.acc[["vista_inactive_all_tissues"]],
         vista_ts= act.acc[["vista_ts"]],
         evo.ini= act.acc[["act_acc_design_ini_seq"]][, id:= .I],
         evo.act= act.acc[["act_acc_design_act_19"]][, id:= .I],
         evo.act.acc= act.acc[["act_acc_design_act_acc_19"]][, id:= .I]),
    ledidi.input,
    ledidi.design
  )
  obj <- rbindlist(obj, idcol = "label", fill = T)
  obj[, width:= nchar(sequence)]
  obj[is.na(selected), selected:= FALSE]
  
  # Import predictions activities for ledidi sequences ----
  pred <- readRDS(
    paste0("/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/Rdata/final_enhancer_selection/all_",
           tiss, "_prediction_score_with_blast.rds"))
  pred <- pred[type=="ledidi_designed_12.0_14.0"]
  pred[, id:= tstrsplit(location, ":", keep= 2)]
  pred[, label:= "ledidi_12_14"]
  pred <- pred[, .(label, id,
                   acc_heart= heart_access, act_heart= heart_sigmoid,
                   acc_limb= limb_access, act_limb= limb_sigmoid,
                   acc_midbrain= midbrain_access, act_midbrain= midbrain_sigmoid,
                   blast= hg19_blast+mm10_blast+vista_blast>0)]
  
  # Add predicted activities for vista test set ----
  meta <- data.table(tissModel= c("heart", "limb", "midbrain"))
  meta <- meta[, .(fold= c("fold01", "fold02", "fold03")), tissModel]
  meta <- meta[, .(rep= c("rep1", "rep2")), .(tissModel, fold)]
  dir <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/result/model/VISTA_model/model1_bulkATAC_tsx3Aug_2xBal_noW/"
  meta <- meta[, {
    folder1 <- paste0(dir, tissModel, "/results_", fold, "_", tissModel, "_DeepSTARR2_", rep, "/sigmoid_prediction/")
    folder2 <- paste0(dir, tissModel, "/results_", fold, "_", tissModel, "_DeepSTARR2_", rep, "/evolutional_enhancer_design/VISTA/Ledidi_two_model_design/predict_access/", tiss, "/")
    .(predicted_act_vista= list.files(folder1, paste0("sequences_test_sigmoid_", tiss, "_predicted_sigmoid_score.txt"), full.names= T),
      sequence_ids= list.files(folder2, paste0("test.fasta_predictions_enhancer_Model.txt"), full.names= T))
  }, .(tissModel, fold, rep)]
  pred2 <- meta[, cbind(fread(predicted_act_vista), fread(sequence_ids)), .(tissModel, predicted_act_vista, sequence_ids)]
  pred2[, id:= tstrsplit(location, ":", keep= 1)]
  setnames(pred2, c("score", "Predictions"), c("act", "acc"))
  pred2 <- dcast(pred2, id~tissModel, value.var = list("act", "acc"), fun.aggregate = list(mean, mean))
  pred2[obj, label:= i.label, on= "id"]
  setnames(pred2, function(x) gsub("_mean", "", x))
  
  # Merge predictions and set cutoffs ----
  pred <- rbind(pred, pred2, fill= T)
  
  # Add predictions to object ----
  obj <- merge(obj,
               pred,
               by= c("label", "id"),
               all.x= T)
  
  # Filtering sequences ----
  if(tiss=="heart") {
    obj[, active:=   label=="ledidi_12_14" & act_heart>6]
    obj[, specific:= label=="ledidi_12_14" & act_midbrain<0 & act_limb<0]
    obj[, selected:= label=="ledidi_12_14" & act_heart>6 & act_midbrain<0 & act_limb<0 & !(blast)]
  }
  if(tiss=="midbrain") {
    obj[, active:=   label=="ledidi_12_14" & act_midbrain>7]
    obj[, specific:= label=="ledidi_12_14" & act_heart<0 & act_limb<0]
    obj[, selected:= label=="ledidi_12_14" & act_midbrain>7 & act_heart<0 & act_limb<0 & !(blast)]
  }
  if(tiss=="limb") {
    obj[, active:=   label=="ledidi_12_14" & act_limb>5]
    obj[, specific:= label=="ledidi_12_14" & act_heart<0 & act_midbrain<0]
    obj[, selected:= label=="ledidi_12_14" & act_limb>5  & act_heart<0 & act_midbrain<0 & !(blast)]
  }
  
  # Save clean object
  saveRDS(obj,
          paste0("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/Rdata/final_designed_enhancer_sequences_", tiss, ".rds"))
}