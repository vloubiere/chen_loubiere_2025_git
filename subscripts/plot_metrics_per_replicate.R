setwd("/groups/stark/vloubiere/projects/DeepATAC_shenzhi/")
require(vlfunctions)

# Import metrics ----
ATAC <- readRDS("db/statistics/model1/ATAC_PCC_replicates.rds")
ATAC <- ATAC[variable=="glob_PCC_tissue" | set=="test"]
ATAC[, set:= switch(set,
                    "chr18"= "Rdm chr18",
                    "test"= "Sel. regions"), set]

VISTA <- readRDS("db/statistics/model1/VISTA_AUC_mPCC_PPV_replicates.rds")
VISTA <- VISTA[variable %in% c("rocAUC", "prAUC", "mcPCC", "predict_cutoff", "PPV_at_cutoff")]
VISTA[, set:= switch(set,
                     "chr18"= "Shared test",
                     "test"= "Fold-spe."), set]

dat <- rbindlist(list(ATAC= ATAC,
                      VISTA= VISTA),
                 idcol = "type")
dat[, set:= factor(set, c("Sel. regions", "Rdm chr18", "Shared test", "Fold-spe."))]
dat[, model:= factor(model, fread("Rdata/model1_ordering.txt", header= F)[[1]])]
setorderv(dat,
          c("type", "set", "variable", "tissue", "model", "fold", "replicate"))


# Compute groups for plotinh ----
Cc <- c("tomato", "cornflowerblue", "limegreen")
dat[, col:= Cc[.GRP], fold]
dat[, line:= .GRP, .(type, set, variable)]
dat[, endLine:= line==max(line)]
dat[, column:= .GRP, tissue]
dat[, endColumn:= column==max(column)]
ncol <-  uniqueN(dat[, .(tissue)])
nrow <-  uniqueN(dat[, .(type, set, variable)]) # Leave space for labels

# Plot ----
dat[, max:= max(unlist(value), na.rm = T), .(type, set, variable)] # compute max
pdf("pdf/model1_performance_metrics_per_replicate.pdf",
    width = 1.6*ncol,
    height = 0.5*nrow)
vl_par(mfrow= c(nrow, ncol),
       mai= c(0, .125, .1, 0),
       omi= c(0.5, 1.25, .25, 0),
       font.main= 1)
# dat[tissue=="forebrain" & variable=="glob_PCC_tissue" & set=="Selected regions", {
dat[, {
  # Barplots
  pl <- .SD[, .(mean= mean(value, na.rm= T), values= .(value), col= .(col)), model]
  bar <- vl_barplot(height= pl$mean,
                    individual.var = pl$values,
                    ind.col = adjustcolor(unlist(pl$col), .5),
                    ylim= c(0, max),
                    ind.cex = .3,
                    border= NA)
  # Legend
  if(column==1)
    text(par("usr")[1],
         mean(par("usr")[3:4]),
         paste0(type, " ", set, "\n", variable),
         pos= 2,
         xpd= NA,
         offset= 2)
  if(line==1)
    title(main= tissue, line= 0.75, xpd= NA)
  if(endLine)
    vl_tilt_xaxis(bar, labels = gsub("^model1_bulkATAC_", "", levels(pl$model)), cex= .4)
}, .(type, set, variable, max, tissue, line, endLine, column, endColumn)]
dev.off()