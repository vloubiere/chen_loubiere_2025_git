setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# For each tissue and dataset ----
output.file <- "db/contributions/mean_contrib_per_motif_instance.rds"
if(!file.exists(output.file)) {
  
  # Metadata
  meta <- data.table(contrib.file= list.files("db/contributions/", recursive = T, full.names = T))
  meta[, dataset:= ifelse(grepl("ATAC", contrib.file), "accessibility", "activity")]
  meta[, tissue:= tstrsplit(basename(contrib.file), "_", keep= 6)]
  meta[, tissue:= gsub(".rds$", "", tissue)]
  
  # Import motif positions
  mot <- readRDS("db/motifs/non_redundant_motif_pos_test_set_paper_subject_0.0001_collapsed.rds")
  
  # Extract contributions per motif
  res <- meta[, {
    
    # Import contribution scores
    contrib <- readRDS(contrib.file)
    contrib <- contrib[strand=="+"]
    contrib[, seq.length:= nchar(seq)]
    if(!"genome" %in% names(contrib))
      contrib[, genome:= "mm10"]
    
    # Compute overlaps
    ov <- overlapBed(mot, contrib)
    ov[, motif:= mot$name[idx.a]]
    ov[, seqnames:= mot$seqnames[idx.a]]
    
    # Make sure the genome is the same
    ov[, genome.a:= mot$genome[idx.a]]
    ov[, genome.b:= contrib$genome[idx.b]]
    ov <- ov[genome.a==genome.b]
    
    # Compute motif start and end position
    ov[, start:= overlap.start-contrib$start[idx.b]+1]
    ov[, end:= start+overlap.width-1]
    ov[, seq.length:= contrib$seq.length[idx.b]]
    
    # Extract mean contrib per motif/control window
    ov[, contrib.mot:= {
      mean(contrib$score[idx.b][[1]][start:end])
    }, .(idx.b, start, end)]
    
    # Sum upstream control
    ov[start>1, upstream.end:= start-1]
    ov[start>1, upstream.start:= ifelse(upstream.end-100>1, upstream.end-100, 1)]
    ov[start>1, upstream.width:= upstream.end-upstream.start+1]
    ov[start>1, upstream.contrib:= {
      sum(contrib$score[idx.b][[1]][upstream.start:upstream.end])
    }, .(idx.b, upstream.start, upstream.end)]
    
    # Sum downstream control
    ov[end<seq.length, downstream.start:= end+1]
    ov[end<seq.length, downstream.end:= ifelse(downstream.start+100>seq.length, seq.length, downstream.start+100)]
    ov[end<seq.length, downstream.width:= downstream.end-downstream.start+1]
    ov[end<seq.length, downstream.contrib:= {
      sum(contrib$score[idx.b][[1]][downstream.start:downstream.end])
    }, .(idx.b, downstream.start, downstream.end)]
    
    # Mean control window
    ov[, contrib.ctl:= (upstream.contrib+downstream.contrib)/(upstream.width+downstream.width)]
    
    # Return
    ov[, .(motif, contrib.mot, contrib.ctl)]
  }, .(tissue, dataset)]
  
  # Save
  saveRDS(res, output.file)
} else
  res <- readRDS(output.file)

