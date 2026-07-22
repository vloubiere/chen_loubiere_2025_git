setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# List observed files ----
dat <- data.table(TL= paste0("TL", c(12, 9, 13)))
dat <- dat[, .(obs.file= list.files(paste0("db/observed/bulkATAC/", TL, "/"), recursive = T, full.names = T)), TL]
dat[, tissue:= tstrsplit(obs.file, "/", keep= 6)]
dat[, c("fold", "set"):= tstrsplit(basename(obs.file), "_|[.]", keep= c(1, 4))]
dat[, c("active", "inactive"):= {
  .c <- fread(obs.file)
  .(sum(.c$score), sum(!.c$score))
}, obs.file]
dat[, total:= active+inactive]

# Save
cols <- c("active", "inactive", "total")
dat[, (cols):= lapply(.SD, formatC, big.mark = ","), .SDcols= cols]
fwrite(dat[, !"obs.file"], "db/statistics/chrom_TL/number_active_sequences.txt", sep= "\t", quote= F, na = NA)
