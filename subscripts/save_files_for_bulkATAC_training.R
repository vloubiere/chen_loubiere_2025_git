setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
source("git_deepATAC/function/augmentation_function_tiling_sliding_window.R")
require(vlfunctions)

# Import models metadata ----
meta <- readxl::read_xlsx("Rdata/metadata_ATACSeq_models.xlsx", skip = 4)
meta <- as.data.table(meta)

# Create output files ----
meta <- meta[, .(tissue= c("heart", "forebrain", "hindbrain", "midbrain", "neuralTube", "limb")), (meta)]
meta <- meta[, .(fold= paste0("fold0", 1:3)), (meta)]
meta <- meta[, .(replicate= paste0("rep", 1:2)), (meta)]
meta <- meta[, .(set= c("training", "validation", "test")), (meta)]
# bed
meta[, bed_VISTA:= paste0("db/bed/", dataset, "/VISTA/", tissue, "/", fold, "_sequences_", set,".bed")]
meta[, bed_ATAC:= paste0("db/bed/", dataset, "/ATAC/", tissue, "/", augmentation, "/", balancing, "/", fold, "_sequences_", set,".bed")]
# fasta
meta[, fa_VISTA:= paste0("db/fasta/", dataset, "/VISTA/", tissue, "/", fold, "_sequences_", set,".fa")]
meta[, fa_ATAC:= paste0("db/fasta/", dataset, "/ATAC/", tissue, "/", augmentation, "/", balancing, "/", fold, "_sequences_", set,".fa")]
# Observed (0/1 labels for VISTA, coverage for ATAC)
meta[, obs_VISTA:= paste0("db/observed/", dataset, "/VISTA/", tissue, "/", fold, "_sequences_activity_", set,".txt")]
meta[, obs_ATAC:= paste0("db/observed/", dataset, "/ATAC/", tissue, "/", augmentation, "/", balancing, "/", fold, "_sequences_activity_", set,".txt")]

# Save processed metadata ----
formula <- as.formula(paste(paste(names(meta)[1:12], collapse = "+"), "~ set"))
meta <- dcast(meta,
              formula,
              value.var = list("bed_VISTA", "bed_ATAC", "fa_VISTA", "fa_ATAC", "obs_VISTA", "obs_ATAC"))
saveRDS(meta,
        "Rdata/metadata_ATACSeq_models_processed.rds")

# Create missing directories ----
dirs <- unique(dirname(unlist(meta[, bed_VISTA_test:obs_ATAC_validation])))
if(any(!dir.exists(dirs)))
  sapply(dirs, function(x) dir.create(x, showWarnings = F, recursive = T))

# Import augmented dataset ----
vista <- readRDS("db/augmentation/bulkATAC/vistaTiles_1001_50.rds")
GO <- readRDS("db/augmentation/bulkATAC/globallyOpen_10x.rds")
TS <- readRDS("db/augmentation/bulkATAC/tissueSpecific_34x.rds")
GC <- readRDS("db/augmentation/bulkATAC/globallyClosed_10x.rds")
GCnoAUG <- readRDS("db/augmentation/bulkATAC/globallyClosed_noAUG.rds")
# Compute 1x(TS+GO+SC after 10x augmentation)
NsampleGC <- length(unique(c(TS$peakID, GO$peakID)))*10

# Import folds ----
folds <- readRDS("db/folds/bulkATAC_folds.rds")

# VISTA tiles ----
overwrite <- F
meta[, {
  if(any(!file.exists(unlist(.BY)[-c(1,2)])) | overwrite)
  {
    # Select VISTA tiles ID
    IDs <- folds[group=="vista", c("peakID", fold), with= F]
    setnames(IDs,
             c("peakID", "set"))
    
    # Merge
    vista[IDs, set:= i.set, on= "peakID"]
    .c <- merge(IDs,
                vista[, c("peakID", "ID", tissue, "seqnames", "start", "end", "strand", "seq"), with= F],
                by= "peakID")
    setnames(.c,
             tissue,
             "score")
    .c[, class:= ifelse(score, "active", "inactive")]
    .c[, ID:= paste0(ID, "__", class)]
    
    # Balance training and validation dataset ----
    .c[, Nact:= sum(score), set]
    .c[, Ntiles:= .N, peakID]
    .c[, prob:= max(Ntiles)/Ntiles]
    sel <- .c[class=="inactive", ID[sample(.N, Nact/5*6, prob= prob)], .(Nact, set)]$V1
    .c <- .c[set %in% c("test", "test.shared") | class=="active" | ID %in% sel]
    
    # Save bed files
    if(!file.exists(bed_VISTA_training) | overwrite)
      rtracklayer::export(GRanges(na.omit(.c[set=="training", .(seqnames, start, end, strand, name= ID, score)])),
                          bed_VISTA_training)
    if(!file.exists(bed_VISTA_validation) | overwrite)
      rtracklayer::export(GRanges(na.omit(.c[set=="validation", .(seqnames, start, end, strand, name= ID, score)])),
                          bed_VISTA_validation)
    if(!file.exists(bed_VISTA_test) | overwrite)
      rtracklayer::export(GRanges(na.omit(.c[grepl("^test", set), .(seqnames, start, end, strand, name= ID, score)])),
                          bed_VISTA_test)
    
    # Save fasta files
    if(!file.exists(fa_VISTA_training) | overwrite)
      seqinr::write.fasta(sequences = as.list(.c[set=="training", seq]), 
                          names = .c[set=="training", ID], 
                          file.out = fa_VISTA_training)
    if(!file.exists(fa_VISTA_validation) | overwrite)
      seqinr::write.fasta(sequences = as.list(.c[set=="validation", seq]), 
                          names = .c[set=="validation", ID], 
                          file.out = fa_VISTA_validation)
    if(!file.exists(fa_VISTA_test) | overwrite)
      seqinr::write.fasta(sequences = as.list(.c[grepl("^test", set), seq]), 
                          names = .c[grepl("^test", set), ID], 
                          file.out = fa_VISTA_test)
    
    # Save observed files
    if(!file.exists(obs_VISTA_training) | overwrite)
      fwrite(.c[set=="training", .(class, ID, score)],
             obs_VISTA_training,
             sep= "\t",
             na = NA)
    if(!file.exists(obs_VISTA_validation) | overwrite)
      fwrite(.c[set=="validation", .(class, ID, score)],
             obs_VISTA_validation,
             sep= "\t",
             na = NA)
    if(!file.exists(obs_VISTA_test) | overwrite)
      fwrite(.c[grepl("^test", set), .(class, ID, score)],
             obs_VISTA_test,
             sep= "\t",
             na = NA)
  }
  print(paste0(.GRP, "/", .NGRP))
}, .(tissue, fold,
     bed_VISTA_training, bed_VISTA_validation, bed_VISTA_test,
     fa_VISTA_training, fa_VISTA_validation, fa_VISTA_test,
     obs_VISTA_training, obs_VISTA_validation, obs_VISTA_test)]

# ATAC peaks ----
meta[, {
  if(any(!file.exists(unlist(.BY)[-c(1,2,3,4)])) | overwrite)
  {
    
    # Select ATAC-seq peaks ID
    IDs <- folds[group=="ATAC", c("peakID", fold), with= F]
    setnames(IDs,
             c("peakID", "set"))
    
    # Globally Open
    GO[IDs, set:= i.set, on= "peakID"]
    .GO <- merge(IDs,
                 GO[, c("peakID", "ID", "shift", paste0(tissue, ".cov"), "seqnames", "start", "end", "strand", "seq", "class", tissue), with= F],
                 by= "peakID")
    setnames(.GO,
             c(tissue, paste0(tissue, ".cov")),
             c("label", "score"))
    .GO[, ID:= paste0(ID, "__", ifelse(label==1, "open", "closed"))]
    if(balancing=="nogo")
      .GO <- .GO[0]
    
    # Tissue Specific
    TS[IDs, set:= i.set, on= "peakID"]
    .TS <- merge(IDs,
                 TS[, c("peakID", "ID", "shift", paste0(tissue, ".cov"), "seqnames", "start", "end", "strand", "seq", "class", tissue), with= F],
                 by= "peakID")
    setnames(.TS,
             c(tissue, paste0(tissue, ".cov")),
             c("label", "score"))
    .TS[, ID:= paste0(ID, "__", ifelse(label==1, "open", "closed"))]
    .TS[, class:= ifelse(label==1, "tissueSpecific", "specificClosed")]
    .TS <- if(augmentation=="tsx3") # TS are all taken (not specificClosed)
      .TS[label==1 | shift %in% seq(-400, 400, 200)] else
        .TS[shift %in% seq(-400, 400, 200)] # Regular 10X
    
    # Globally Closed
    IDs <- folds[group=="ctl", c("peakID", fold), with= F]
    setnames(IDs,
             c("peakID", "set"))
    .GC <- if(augmentation=="nogc")
    {
      GCnoAUG[IDs, set:= i.set, on= "peakID"]
      merge(IDs,
            GCnoAUG[, c("peakID", "ID", "shift", paste0(tissue, ".cov"), "seqnames", "start", "end", "strand", "seq", "class", tissue), with= F],
            by= "peakID")
    }else
    {
      GC[IDs, set:= i.set, on= "peakID"]
      merge(IDs,
            GC[, c("peakID", "ID", "shift", paste0(tissue, ".cov"), "seqnames", "start", "end", "strand", "seq", "class", tissue), with= F],
            by= "peakID")
    }
    setnames(.GC,
             c(tissue, paste0(tissue, ".cov")),
             c("label", "score"))
    .GC[, ID:= paste0(ID, "__", ifelse(label==1, "open", "closed"))]
    
    # Balancing
    .GC <- if(balancing=="1x")
    {
      set.seed(1)
      .GC[sample(.N, NsampleGC)]
    }else
    {
      set.seed(1)
      .GC[sample(.N, NsampleGC*2)]
    }
    
    # Combine and print stats
    cmb <- rbind(.GO, .TS, .GC)
    fwrite(cmb[, .N, class],
           paste0(dirname(obs_ATAC_training), "/classes_stats.txt"),
           sep= "\t")
    
    # Save bed files
    if(!file.exists(bed_ATAC_training) | overwrite)
      rtracklayer::export(GRanges(na.omit(cmb[set=="training", .(seqnames, start, end, strand, name= ID, score)])),
                          bed_ATAC_training)
    if(!file.exists(bed_ATAC_validation) | overwrite)
      rtracklayer::export(GRanges(na.omit(cmb[set=="validation", .(seqnames, start, end, strand, name= ID, score)])),
                          bed_ATAC_validation)
    if(!file.exists(bed_ATAC_test) | overwrite)
      rtracklayer::export(GRanges(na.omit(cmb[grepl("^test", set), .(seqnames, start, end, strand, name= ID, score)])),
                          bed_ATAC_test)
    
    # Save fasta files
    if(!file.exists(fa_ATAC_training) | overwrite)
      seqinr::write.fasta(sequences = as.list(cmb[set=="training", seq]), 
                          names = cmb[set=="training", ID], 
                          file.out = fa_ATAC_training)
    if(!file.exists(fa_ATAC_validation) | overwrite)
      seqinr::write.fasta(sequences = as.list(cmb[set=="validation", seq]), 
                          names = cmb[set=="validation", ID], 
                          file.out = fa_ATAC_validation)
    if(!file.exists(fa_ATAC_test) | overwrite)
      seqinr::write.fasta(sequences = as.list(cmb[grepl("^test", set), seq]), 
                          names = cmb[grepl("^test", set), ID], 
                          file.out = fa_ATAC_test)
    
    # Save observed files
    if(!file.exists(obs_ATAC_training) | overwrite)
      fwrite(cmb[set=="training", .(class, ID, score)],
             obs_ATAC_training,
             sep= "\t",
             na = NA)
    if(!file.exists(obs_ATAC_validation) | overwrite)
      fwrite(cmb[set=="validation", .(class, ID, score)],
             obs_ATAC_validation,
             sep= "\t",
             na = NA)
    if(!file.exists(obs_ATAC_test) | overwrite)
      fwrite(cmb[grepl("^test", set), .(class, ID, score)],
             obs_ATAC_test,
             sep= "\t",
             na = NA)
  }
  print(paste0(.GRP, "/", .NGRP))
}, .(augmentation, balancing, tissue, fold, 
     bed_ATAC_training, bed_ATAC_validation, bed_ATAC_test,
     fa_ATAC_training, fa_ATAC_validation, fa_ATAC_test,
     obs_ATAC_training, obs_ATAC_validation, obs_ATAC_test)]
