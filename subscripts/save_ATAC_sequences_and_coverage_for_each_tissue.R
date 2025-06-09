setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import ATAC-Seq folds ----
ATAC <- readRDS("Rdata/folds_training_valid_test_sets_cov_act.rds")[modelType=="ATAC"]

# For each tissue, retrieve the accessibility of each class  ----
tiss <- data.table(tissue= c("heart", "limb", "forebrain", "hindbrain", "midbrain", "neuralTube"))
tiss <- tiss[, .(class= unique(ATAC$class)), tissue]
tiss[, label:= "closed"]
tiss[class=="globallyOpen", label:= "open"]
tiss[tissue=="heart" & grepl("heart", class, ignore.case = T), label:= "open"]
tiss[tissue=="limb" & grepl("limb", class, ignore.case = T), label:= "open"]
tiss[tissue=="forebrain" & grepl("forebrain|panBrain", class, ignore.case = T), label:= "open"]
tiss[tissue=="hindbrain" & grepl("hindbrain|panBrain", class, ignore.case = T), label:= "open"]
tiss[tissue=="midbrain" & grepl("midbrain|panBrain", class, ignore.case = T), label:= "open"]
tiss[tissue=="neuralTube" & grepl("neural|panBrain", class, ignore.case = T), label:= "open"]

# Combined folds and tissues ----
folds <- melt(ATAC,
              id.vars = c("ID", "sequence", "class"),
              measure.vars = patterns("^fold"),
              variable.name = "fold")
folds <- merge(folds,
               tiss,
               by= "class",
               allow.cartesian= T)

# Subsample globally closed to match opened regions ----
folds[, Nopen:= sum(label=="open"), .(fold, tissue, value)]
folds <- folds[, {
  set.seed(.GRP)
  rbind(.SD[class!="globallyClosed"],
        .SD[class=="globallyClosed"][sample(.N, Nopen*2)])
}, .(fold, tissue, value, Nopen)]

# Simplify classes and check ----
folds[class!="globallyClosed" & label=="closed", class:= "specificClosed"]
check <- folds[, .N, .(tissue, fold, class)]
checkT <- folds[, .N, .(tissue, fold, value)]

# Retrieve coverage in the corresponding tissue ----
cov <- melt(ATAC[, .(ID, heart, limb, forebrain, hindbrain, midbrain, neuralTube)],
            id.vars = "ID")
folds[cov, score:= i.value, on= c("ID", "tissue==variable")]
folds[ATAC, seq:= i.sequence, on= "ID"]

# Save fasta files ----
folds[, fa_file:= {
  .dir <- paste0("db/fasta/ATAC/", tissue, "/")
  dir.create(.dir, showWarnings = F, recursive = T)
  paste0(.dir, fold, "_sequences_", value, ".fa")
}, .(tissue, fold, value)]
folds[, {
  seqinr::write.fasta(as.list(seq),
                      ID,
                      fa_file)
  print(fa_file)
}, fa_file]

# Save txt files ----
folds[, txt_file:= {
  .dir <- paste0("db/scores/ATAC/", tissue, "/")
  dir.create(.dir, showWarnings = F, recursive = T)
  paste0(.dir, fold, "_sequences_activity_", value, ".txt")
}, .(tissue, fold, value)]
folds[, {
  fwrite(.SD,
         txt_file,
         col.names = T,
         sep= "\t")
  print(txt_file)
}, txt_file, .SDcols= c("class", "ID", "score")]