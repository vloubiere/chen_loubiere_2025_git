setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[dataset=="activity" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set=="test"]

# Output files ----
meta[, mean_contrib_file:= file.path(
  "db/contributions/VISTA/",
  paste0(dataset, "_model_mean_contrib_scores_", tissue, ".rds")
), .(dataset, tissue)]

# Compute rolling mean ----
overwrite <- FALSE
meta[, {
  if(overwrite | !file.exists(mean_contrib_file)) {
    # Import data
    contrib <- .SD[, importContrib(contrib_file, bed = bed_file, fa = fa_file, FUN = colSums), .(contrib_file, bed_file, fa_file)]
    contrib[, genome:= ifelse(grepl("^hs", name), "hg38", "mm10")]
    # For each genome and strand
    res <- contrib[, {
      # Collapsed window indexes
      idx <- collapseBed(.SD, return.idx.only = T)
      # Melt data
      .SD[, {
        start <- seq(coor, coor+1000)
        .(start= if(strand=="+") start else rev(start),
          base= strsplit(seq, "")[[1]],
          score= unlist(contrib.score))
      }, .(idx, seqnames, coor= start)]
    }, .(contrib_file, genome, strand)]
    # Mean score per position
    res <- res[, .(score= mean(score)), keyby= .(genome, idx, seqnames, start, base, strand)]
    # Collapse
    res <- res[, .(start= start[1], seq= paste0(base, collapse= ""), score= list(score)), keyby= .(genome, idx, seqnames, strand)]
    res[, end:= start+nchar(seq)-1]
    # Simplify
    res <- res[, .(genome, seqnames, start, end, strand, seq, score)]
    # Collapse
    saveRDS(res, mean_contrib_file)
  }
  print(paste(tissue, "DONE"))
}, .(tissue, mean_contrib_file)]

