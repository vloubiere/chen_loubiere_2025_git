setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Import metadata
meta <- readRDS("Rdata/paper_metadata_v1.rds")
meta <- meta[dataset=="accessibility" & ID=="model1_bulkATAC_tsx3Aug_2xBal_noW" & set=="chr18"]

# Create accessibility predicted bedgraph files ----
meta[, bedgraph.file:= {
  # Output file path
  bedgraph.file <- paste0("db/bedgraph/predicted/", tissue, "_predicted_accessibility_chr18.bedgraph")
  if(!file.exists(bedgraph.file)) {
    # Import all folds/replicates predictions
    .c <- .SD[, {
      fread(pred_file)
    }, .(fold, replicate, pred_file)]
    # Only keep tiles on the + strand (- is redundant)
    .c <- .c[grepl("_\\+$", location)]
    .c[, location:= gsub("_\\+$", "", location)]
    # Compute mean across all fold/replicates
    .c <- .c[, .(Predictions= mean(Predictions)), .(location)]
    # Extract bed coordinates
    .c <- cbind(importBed(.c$location), .c[, .(Predictions)])
    # Mean across overlapping bins
    .c[, end:= end-1]
    bins <- binBed(.c, bins.width = 200)
    # Unlog signal (to compare with ATAC-Seq)
    mean <- bins[, .(score= 2^mean(Predictions)), .(seqnames, start, end)]
    # Save bedraph file
    rtracklayer::export(mean, bedgraph.file)
  }
  bedgraph.file
}, tissue]

# Bedgraph to bigwig ----
meta[, {
  # Make command
  cmd <- vlite::cmd_bedgraphToBigwig(
    bdg = bedgraph.file,
    genome = "mm10",
    bw.output.folder = "db/bw/predicted/"
  )
  # Submit
  vl_submit(cmd,
            overwrite = FALSE)
}, bedgraph.file]