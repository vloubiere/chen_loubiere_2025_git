setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")

# Make metadata ----
meta <- data.table(bed_file= list.files("db/bed/bulkATAC/", "_test.bed$", recursive = T, full.names = T))
meta[, dataset:= tstrsplit(bed_file, "/", keep= 5)]
meta[dataset=="ATAC", c("tissue", "augmentation", "balancing"):= tstrsplit(bed_file, "/", keep= 6:8)]
meta[dataset=="VISTA", tissue:= tstrsplit(bed_file, "/", keep= 6)]
meta[, fold:= tstrsplit(basename(bed_file), "_", keep= 1)]

# Test set is already unique for VISTA (no balancing, only one) ----
meta <- meta[dataset!="VISTA"]

# Problem: in our old design, the test set differed between different strategies -> PCC can't be compared ----
# So for consistency, I ill now just copy and paste the test set of our favorite design in all folders
meta[, test_bed_file_best:= bed_file[augmentation=="tsx3" & balancing=="2x"], .(tissue, fold)]
meta[, copy_name:= gsub("_test.bed$", "_testBestDesign.bed", bed_file)]
meta[, file.copy(test_bed_file_best, copy_name, overwrite = T), .(test_bed_file_best, copy_name)]

# Before -> 5 different test sets per fold
meta[, {
  .c <- lapply(unique(bed_file), fread, sel= 1:3)
  cache <- sapply(.c, digest)
  .(N_different_test_sets= length(unique(cache)))
}, .(dataset, fold, tissue)]
# After copy -> 1 best test per fold
meta[, {
  .c <- lapply(unique(copy_name), fread, sel= 1:3)
  cache <- sapply(.c, digest)
  .(N_different_test_sets= length(unique(cache)))
}, .(dataset, fold, tissue)]

# Sanity check no overlaps ----
meta[, {
  bestTest <- importBed(copy_name)
  valid <- importBed(gsub("_testBestDesign.bed$", "_validation.bed", copy_name))
  train <- importBed(gsub("_testBestDesign.bed$", "_training.bed", copy_name))
  stopifnot(intersectBed(bestTest, valid)==0)
  stopifnot(intersectBed(bestTest, train)==0)
  stopifnot(intersectBed(valid, train)==0)
  print(paste(tissue, fold, "OK"))
}, copy_name]
