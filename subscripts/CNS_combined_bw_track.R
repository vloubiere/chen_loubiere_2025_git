setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

cmd <- vlite::cmd_mergeBigwig(
  bw = c(
    "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/clean_version/db/bw/forebrain_merged.CPM.bigwig",
    "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/clean_version/db/bw/hindbrain_merged.CPM.bigwig",
    "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/clean_version/db/bw/midbrain_merged.CPM.bigwig",
    "/groups/stark/shenzhi.chen/projects/transferLearningMammalianEnhancerDesign202408/clean_version/db/bw/neural-tube_merged.CPM.bigwig"
  ),
  output.prefix = "CNS_merged_CPM",
  bw.output.folder = "db/bw/observed/",
  genome = "mm10"
)

vl_submit(cmd, mem = 128) # My method imports the bedgraph into R -> lots of mem needed here :()

