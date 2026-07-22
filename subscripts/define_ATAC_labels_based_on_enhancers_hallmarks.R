setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite-dev/")

# Import metadata ----
meta <- readRDS("Rdata/paper_metadata_v3.rds")
meta <- meta[set %in% c("test", "validation", "training")]
meta <- meta[dataset %in% "accessibility"]
meta <- meta[tissue %in% c("heart", "midbrain", "limb")]
meta <- meta[ID=="model1_bulkATAC_tsx3Aug_2xBal_noW"]

# ChIPSeq peaks ----
peaks <- readRDS("Rdata/revision_ChIPseq_peaks.rds")
peaks <- peaks[signalValue>5]

# ATAC peaks ----
ATAC <- readRDS("db/peaks/bulkENCODE_confident_ATAC_peaks.rds")
# Tissue-specific
ATAC <- readRDS("db/peaks/bulkENCODE_confident_ATAC_peaks.rds")
ATAC <- melt(ATAC, c("seqnames", "start", "end"), c("limb", "midbrain", "heart"))
ATAC <- ATAC[value==1]

# Promoters ----
if(!exists("ext_prom")) {
  gtf <- rtracklayer::import("../../genomes/Mus_musculus/GENCODE/gencode.vM25.annotation.gtf.gz")
  gtf <- as.data.table(gtf)
  gtf <- gtf[type=="transcript"]
  ext_prom <- resizeBed(gtf[, .(seqnames, start, end)], "start", 10000, 10000)
  ext_prom <- unique(ext_prom)
}

# Tests ----
meta[, {
  
  # Import bed file ----
  .c <- importBed(bed_file)
  
  # Select only center tiles (no augmentation)
  .c[, c("ID", "tile", "label"):= tstrsplit(name, "_|__")]
  .c <- .c[tile==0]
  
  # Compute overlaps with enhancer marks ----
  .c$K27Ac_ov <- covBed(.c, peaks[.BY, on= "tissue"][cdition=="H3K27ac"])>0
  .c$K4me1_ov <- covBed(.c, peaks[.BY, on= "tissue"][cdition=="H3K4me1"])>0
  .c$K4me3_ov <- covBed(.c, peaks[.BY, on= "tissue"][cdition=="H3K4me3"])>0
  
  # Compute overlaps with genomic annotations ----
  .c$ext_prom <- covBed(.c, ext_prom)>0
  .c$CTCF_ov <- covBed(.c, peaks[cdition=="CTCF"])>0 # Union of all tissues
  
  # Tissue-specific ATAC-seq peaks ----
  .c$ATAC_ov <- covBed(.c, ATAC[variable==tissue])>0
  .c$ATAC_other_tissue <- covBed(.c, ATAC[variable!=tissue])>0
  
  # Import observed file and select center tiles ----
  obs <- fread(obs_file)
  obs <- obs[ID %in% .c$name]
  
  # Compute TL enhancer labels ----
  obs[, ACC:= ID %in% .c[(ATAC_ov), name]] # Accessible in the tissue
  obs[, TL8:= ACC & ID %in% .c[(K27Ac_ov), name]]
  obs[, TL9:= ACC & ID %in% .c[(K27Ac_ov & K4me1_ov), name]]
  obs[, TL10:= ACC & ID %in% .c[(K27Ac_ov & K4me1_ov & !K4me3_ov), name]]
  obs[, TL11:= ACC & ID %in% .c[(!ext_prom), name]]
  obs[, TL12:= ACC & ID %in% .c[(!ext_prom) & (!CTCF_ov), name]]
  obs[, TL13:= ACC & ID %in% .c[(!ATAC_other_tissue), name]]
  obs[, TL14:= TL10 & TL13] # Combine HTMs and tissue-specific ATAC
  
  # Reformat observed file ----
  obs$ACC <- NULL # Remove accessibility column
  obs <- melt(obs, id.vars= c("class", "ID", "score"))
  obs[, score:= ifelse(value, 1, 0)]
  obs[, class:= ifelse(value, "active", "inactive")]
  obs[, ID.new:= tstrsplit(ID, "__", keep= 1)]
  obs[, ID.new:= paste0(ID.new, "__", class)]
  
  # Compute fraction of positive sequences ----
  obs[, Npos:= sum(score), variable]
  obs[, perc.pos:= Npos/.N*100, variable]
  
  # Down-sample if necessary (~13% of positive in the original VISTA models) ----
  obs[, Nneg.to.sample := round(Npos * (1 - 0.13) / 0.13), by = variable]
  obs[, row.idx:= .I]
  obs[score==0 & grepl("__closed", ID), subsample:= {
    if(perc.pos < 13)
      row.idx %in% sample(row.idx, size = Nneg.to.sample) else
        TRUE
  }, .(perc.pos, variable, Nneg.to.sample)]
  obs <- obs[(score==1) | subsample]
  # Sanity check down-sampling
  obs[, Npos:= sum(score), variable]
  obs[, perc.sub:= Npos/.N*100, variable]
  
  # Import fasta file and select center tiles  ----
  fas <- seqinr::read.fasta(fa_file, as.string = T)
  fas <- fas[names(fas) %in% .c$name] # Select center tiles
  
  # Save observed and fasta files ----
  obs[, {
    
    # Create obs directory
    dir.obs <- file.path("db/observed/bulkATAC", variable, tissue)
    dir.create(dir.obs, showWarnings = F, recursive= T)
    
    # Save observed file
    fwrite(
      .SD[, .(class, ID= ID.new, score)],
      file.path(dir.obs, basename(obs_file)),
      sep = "\t",
      quote = FALSE,
      na = NA
    )
    
    # Create fasta directory
    dir.fasta <- file.path("db/fasta/bulkATAC", variable, tissue)
    dir.create(dir.fasta, showWarnings = F, recursive= T)
    
    # Subsample fasta file
    .f <- fas[ID]
    
    # Save modified fasta file
    seqinr::write.fasta(
      sequences = .f,
      names = ID.new,
      file.out = file.path(dir.fasta, basename(fa_file))
    )
  }, variable]
  
  # Progression
  print(paste(unlist(.BY[1:3])))
}, .(set, tissue, fold, bed_file, obs_file, fa_file)]
