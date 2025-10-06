setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")
require(orthogene)

# Import cluster metadata ----
jeff.meta <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/db/motif/motif_annotations.xlsx"
meta <- readxl::read_xlsx(jeff.meta, 1)
meta <- as.data.table(meta)
meta1 <- readxl::read_xlsx(jeff.meta, 2)
meta1 <- as.data.table(meta1)
meta <- merge(
  x = meta[, .(Cluster_ID, Name)],
  y = meta1[, .(Cluster_ID, Motif)],
  by= "Cluster_ID"
)
# Add cluster to ID
meta[, ID:= paste0(Name, "__", Motif)]
meta[, TF:= tstrsplit(Motif, "_", keep= 1)]
meta[, TF:= lapply(TF, function(x) strsplit(x, "\\+")[[1]])]
meta[, TF:= lapply(TF, function(x) gsub(".mouse$", "", x))]

# Add gene IDs
IDs <- meta[, .(TF= unlist(TF)), ID]
mm <- as.data.table(orthogene::map_genes(IDs$TF, species = "mouse"))
mm <- unique(mm[, .(input, mouse_id= target)])
hs <- as.data.table(orthogene::map_genes(IDs$TF, species = "human"))
hs <- unique(hs[, .(input, human_id= target)])
cmb <- merge(mm, hs)
IDs <- merge(IDs, cmb, by.x= "TF", by.y= "input")[, lapply(.SD, list), ID]
meta[IDs, mouse_id:= i.mouse_id, on= "ID"]
meta[IDs, human_id:= i.human_id, on= "ID"]

# Import 2,179 motifs from Shenzhi ----
pwm.file <- "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/clean_version/db/motif/jeff_non-redundant_pwms_jaspar_all.rds"
mot <- readRDS(pwm.file)

# Subset only seed motifs and rename accordingly ----
sel <- match(meta$Motif, sapply(mot, name))
my.mot.collection <- mot[sel]
for(i in seq_along(my.mot.collection))
  my.mot.collection[[i]]@name <- meta$ID[i]
name(my.mot.collection)

# Save ----
final <- list(meta= meta, pwms_log_odds= my.mot.collection)
saveRDS(final, "/groups/stark/vloubiere/motifs_db/non_redudant_mammals_Jeff_motifs_full.rds")
