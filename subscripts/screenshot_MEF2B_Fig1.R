setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
devtools::load_all("/groups/stark/vloubiere/vlite/")

# Screenshot
pdf("pdf/0_paper/ATAC_screenshot_MEF2B_locus.pdf", height = 4)
vlite::bwScreenshot(bed = "chr8:70,145,000-70,175,000",
                    tracks = "/groups/stark/vloubiere/projects/DeepATAC_shenzhi/db/bw/observed/heart_treat_pileup.bigwig",
                    col = "black",
                    genome = "mm10",
                    ngenes = 5)
dev.off()