setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import VISTA sequences ----
vista <- readRDS("Rdata/folds_training_valid_test_sets_cov_act.rds")[modelType=="VISTA"]

# melt data ----
folds <- melt(vista,
              id.vars = c("ID", "sequence", "class"),
              measure.vars = patterns("^fold"),
              variable.name = "fold")
tiss <- c("heart", "limb", "forebrain", "hindbrain", "midbrain", "neuralTube")
folds <- merge(folds,
               melt(vista,
                    id.vars = "ID",
                    measure.vars = tiss,
                    variable.name = "tissue",
                    value.name = "label"),
               by= "ID",
               allow.cartesian= T)
folds[, label:= as.character(ifelse(label, "Active", "Inactive"))]
folds[, VistaID:= tstrsplit(ID, "[.]", keep= 1)]

# Remove inactive VISTA enhancers that have less then ten tiles (after removing duplicated sequences...)
# remove <- folds[label=="Inactive", .N<10, .(fold, tissue, VistaID)]
# remove <- unique(remove[(V1), VistaID])
# tot <- nrow(folds)
# folds <- folds[!(VistaID %in% remove)]
# print(paste(length(remove), "inactive VISTA enhancers had <10 tiles and were removed, corresponding to", tot-nrow(folds), "tiles"))

# # Subsample inactive VISTA enhancers to 10 tiles each ----
# folds[label=="Inactive" & value %in% c("training", "validation"), select:= {
#   .c <- rep(F, .N)
#   set.seed(.GRP)
#   .c[sample(.N, 10)] <- T
#   .c
# }, .(fold, tissue, VistaID)]
# folds[is.na(select), select:= TRUE]
# tot <- nrow(folds)
# folds <- folds[(select)]
# print(paste(tot-nrow(folds), "inactive VISTA tiles removed"))

# Look at the ratio between active and inactive for each fold/tissue ---
folds[, Nact:= sum(label=="Active"), .(fold, tissue, value)]
folds[, Ninact:= sum(label=="Inactive"), .(fold, tissue, value)]
barplot(folds[value=="training", Ninact/(Nact+Ninact), .(fold, tissue, Nact, Ninact)]$V1,
        ylab= "Fraction of active sequences in training set")
barplot(folds[value=="validation", Ninact/(Nact+Ninact), .(fold, tissue, Nact, Ninact)]$V1,
        ylab= "Fraction of active sequences in training set")

# Save fasta files ----
folds[, nameID:= paste0(ID, "_", label)]
folds[, fa_file:= {
  .dir <- paste0("db/fasta/VISTA/", tissue, "/")
  dir.create(.dir, showWarnings = F, recursive = T)
  paste0(.dir, fold, "_sequences_", value, ".fa")
}, .(tissue, fold, value)]
folds[, {
  seqinr::write.fasta(as.list(sequence),
                      nameID,
                      fa_file)
  print(fa_file)
}, fa_file]

# Save txt files ----
folds[, txt_file:= {
  .dir <- paste0("db/scores/VISTA/", tissue, "/")
  dir.create(.dir, showWarnings = F, recursive = T)
  paste0(.dir, fold, "_sequences_activity_", value, ".txt")
}, .(tissue, fold, value)]
folds[, {
  fwrite(.SD,
         txt_file,
         col.names = T,
         sep= "\t")
  print(txt_file)
}, txt_file, .SDcols= c("class", "nameID", "label")]
